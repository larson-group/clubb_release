"""Restore Native physical continuation groups for every science call; keep tokens exact."""
from pathlib import Path
import re,json,sys,hashlib,textwrap
D=Path(__file__).resolve().parent;R=D.parents[2];sys.path.insert(0,str(R/'clubb_kokkos/src/f2k'))
from generate_derived_types import split_top
from fortran_procedures import find_procedure
O=json.loads((D/'SOURCE_OUTLINE.json').read_text());S=(R/O['source']).read_text();names=O['native_inventory_and_order'];contracts={p['name']:p['arguments'] for p in O['routines']}
contracts['lin_interpolate_two_points']=[vars(a) for a in find_procedure('lin_interpolate_two_points').arguments]
def clean(s):return re.sub(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"',lambda m:re.sub(r'[^\n]',' ',m[0]) if m[0].startswith('/') else m[0],s,flags=re.S)
def toks(s):return re.findall(r'"(?:\\.|[^"\\])*"|\w+|[^\s]',clean(s))
def bal(s,op,l='(',r=')'):
 end=op+1;depth=1
 while depth:depth+=(s[end]==l)-(s[end]==r);end+=1
 return s[op+1:end-1],end
def native_groups(lines,start,op):
 rows=[];depth=1
 for ii,line in enumerate(lines[start:]):
  code=line.partition('!')[0]
  if ii==0:code=code[op+1:]
  chars=[]
  for c in code:
   depth+=(c=='(')-(c==')')
   if depth==0:break
   chars.append(c)
  row=''.join(chars).strip('& \t').rstrip(',').strip()
  if row:rows.append(split_top(row))
  if depth==0:break
 return rows
blocks={p['name']:'\n'.join(S.splitlines()[p['source_lines'][0]-1:p['source_lines'][1]]) for p in O['routines']}
calls={}
for name,block in blocks.items():
 lines=block.splitlines();found=[]
 for i,line in enumerate(lines):
  code=line.partition('!')[0]
  if re.search(r'\b(?:subroutine|function)\b',code,re.I):continue
  for m in re.finditer(r'\b('+ '|'.join(contracts)+r')\s*\(',code,re.I):
   group=native_groups(lines,i,m.end()-1);count=sum(map(len,group))
   minimum=sum(not a['optional'] for a in contracts[m[1]])
   if minimum<=count<=len(contracts[m[1]]):found.append((m[1],group))
 calls[name]=found
report={}
for ext in ('cpp','hpp'):
 p=D/f'remapping_module.{ext}';original=p.read_text();C=clean(original);before=original;edits=[];bodies=[];grouped=0
 if ext=='cpp':
  for m in re.finditer(r'\b('+ '|'.join(names)+r')\s*\(',C):
   a,end=bal(C,m.end()-1);pos=end
   while C[pos].isspace():pos+=1
   if C[pos]!='{':continue
   actual=split_top(a)
   if not all(re.search(r'\s+[*&]?\w+(?:\s*=.*)?$',x) for x in actual):continue
   body,stop=bal(C,pos,'{','}')
   if m[1]=='kmppm' and len(actual)==3:continue
   bodies.append((m[1],pos+1,stop-1,calls[m[1]]))
 else:
  # Thin generic View wrapper follows the Native public signature groups.
  m=re.search(r'\bremap_vals_to_target\s*\(',C);a,end=bal(C,m.end()-1)
  # first occurrence is declaration; second is generic definition.
  m=re.search(r'\bremap_vals_to_target\s*\(',C[end:]);op=end+m.end()-1
  _,finish=bal(C,op);pos=finish
  while C[pos].isspace():pos+=1
  body,stop=bal(C,pos,'{','}')
  raw=blocks['remap_vals_to_target'].splitlines();first=next(i for i,l in enumerate(raw) if 'function remap_vals_to_target(' in l)
  group=native_groups(raw,first,raw[first].index('('));bodies=[('generic',pos+1,stop-1,[('remap_vals_to_target',group)])]
 for owner,start,stop,expected in bodies:
  callees={n.lower() for n,g in expected};actual=[]
  if not callees:continue
  for m in re.finditer(r'\b('+ '|'.join(callees)+r')\s*\(',C[start:stop],re.I):
   op=start+m.end()-1;raw,end=bal(C,op);actual.append((m[1],op,end,split_top(raw)))
  assert [n.lower() for n,g in expected]==[n.lower() for n,o,e,a in actual],(owner,expected,actual)
  for (name,groups),(actual_name,op,end,args) in zip(expected,actual):
   sizes=[len(g) for g in groups];n=sum(sizes)
   tags=[{'in':'In','inout':'InOut','out':'Out','':'In',None:'In'}[a['intent']] for a in contracts[name]][:n]
   if name=='kmppm' and len(args)==3:
    # Fused cell limiter receives scalar dm, one coefficient subview and lmt;
    # the Native itot argument belongs only to the array dispatch boundary.
    flat=[x for row in groups for x in row];remove=2;offset=0;new_sizes=[]
    for size in sizes:
     new_sizes.append(size-(offset<=remove<offset+size));offset+=size
    sizes=[x for x in new_sizes if x];tags.pop(remove);n-=1
   assert len(args)==n,(owner,name,n,args)
   args=[re.sub(r'\s+',' ',a).strip() for a in args]
   line=C.rfind('\n',0,op)+1;indent=re.match(r'\s*',C[line:op])[0]+'    ';out=[];index=0
   for size in sizes:
    row=[];tag=None
    for _ in range(size):
     value=args[index];nexttag=tags[index];index+=1
     if row and (tag!=nexttag or len(indent+', '.join(row+[value])+', // '+nexttag)>100):out.append((row,tag));row=[]
     row.append(value);tag=nexttag
    if row:out.append((row,tag))
   rendered='\n'+'\n'.join(indent+', '.join(row)+(',' if i<len(out)-1 else '')+' // '+tag for i,(row,tag) in enumerate(out))+'\n'+indent[:-4]
   edits.append((op+1,end-1,rendered));grouped+=1
 for a,b,value in sorted(edits,reverse=True):original=original[:a]+value+original[b:]
 assert toks(original)==toks(before),ext
 p.write_text(original)
 report[p.name]={'before_sha256':hashlib.sha256(before.encode()).hexdigest(),'after_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'science_calls_restored_to_Native_groups':grouped,'exact_executable_tokens_unchanged':True,'token_count':len(toks(original))}
(D/'review_call_format_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
