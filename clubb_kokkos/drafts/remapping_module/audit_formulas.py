"""Check PPM assignment trees and every explicit array coordinate, without physics runs."""
from pathlib import Path
import re,sys,json,ast,hashlib
D=Path(__file__).resolve().parent;R=D.parents[2];sys.path.insert(0,str(R/'clubb_kokkos/src/f2k'))
from fortran_procedures import statements,declaration
from generate_derived_types import split_top
S=(R/'src/CLUBB_core/remapping_module.F90').read_text();C=(D/'remapping_module.cpp').read_text()
# Keep strings intact; blank comments without shifting source offsets.
def clean(s):return re.sub(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"',lambda m:re.sub(r'[^\n]',' ',m[0]) if m[0].startswith('/') else m[0],s,flags=re.S)
def balanced(s,op,l='(',r=')'):
 end=op+1;depth=1
 while depth:depth+=(s[end]==l)-(s[end]==r);end+=1
 return s[op+1:end-1],end
def poly(s):
 def add(a,b,scale=1):
  d=dict(a)
  for k,v in b.items():d[k]=d.get(k,0)+v*scale
  return {k:v for k,v in d.items() if v}
 def n(x):
  if isinstance(x,ast.Constant):return {'1':x.value}
  if isinstance(x,ast.Name):return {x.id:1}
  if isinstance(x,ast.UnaryOp):return {k:v*(-1 if isinstance(x.op,ast.USub) else 1) for k,v in n(x.operand).items()}
  if isinstance(x,ast.BinOp):
   a,b=n(x.left),n(x.right)
   if isinstance(x.op,ast.Add):return add(a,b)
   if isinstance(x.op,ast.Sub):return add(a,b,-1)
   if isinstance(x.op,ast.Mult):
    if set(a)<={'1'}:return {k:a.get('1',0)*v for k,v in b.items()}
    if set(b)<={'1'}:return {k:b.get('1',0)*v for k,v in a.items()}
  raise ValueError(s)
 return tuple(sorted(n(ast.parse(s,mode='eval').body).items()))
O=json.loads((D/'SOURCE_OUTLINE.json').read_text());NB={p['name']:'\n'.join(S.splitlines()[p['source_lines'][0]-1:p['source_lines'][1]]) for p in O['routines']}
CB={};CC=clean(C)
for n in ['kmppm','steepz','ppm2m','map1_ppm']:
 for m in re.finditer(r'\b'+n+r'\s*\(',CC):
  args,end=balanced(CC,m.end()-1);pos=end
  while CC[pos].isspace():pos+=1
  if CC[pos]!='{':continue
  body,stop=balanced(CC,pos,'{','}')
  if n=='kmppm' and 'real_t dm' not in args:continue
  CB[n]=body;break
 assert n in CB,n
reports={}
for name in CB:
 meta={}
 for st in statements(NB[name]):
  if '::' in st:meta.update(declaration(st))
 arrays={k for k,v in meta.items() if v.dims}
 # kmppm's old assumed-size declaration is understood by the existing parser.
 assert {'a4','dm'}<=arrays if name=='kmppm' else True
 def expr(s,cpp=False,key=False):
  s=re.sub(r'\s+',' ',s).lower().replace('kokkos::','').replace('_core_rknd','').replace('copysign(','sign(')
  if cpp:
   s=re.sub(r'\bcolumn\b','(native_i-i1)',s)
   s=re.sub(r'\bi\b','(native_i-1)' if name=='map1_ppm' else '(native_i-i1)',s)
  else:s=re.sub(r'\bi\b','native_i',s)
  refs=[]
  pattern=r'\b('+'|'.join(sorted(arrays,key=len,reverse=True))+r')\s*\('
  while m:=re.search(pattern,s):
   raw,end=balanced(s,m.end()-1);idx=split_top(raw)
   if name=='kmppm':
    if m[1]=='dm':replacement='dm';s=s[:m.start()]+replacement+s[end:];continue
    idx=idx[:1]
    if cpp:idx=['('+x+')+1' for x in idx]
   elif cpp:
    if m[1] in {'a4','q4'}:
     idx=['('+idx[2]+')+1','('+idx[0]+')+i1','('+idx[1]+')+1']
    elif name=='map1_ppm' and m[1] in {'pe1','pe2','q1','q2'}:
     idx=['('+x+')+1' for x in idx]
     if m[1] in {'q1','q2'}:idx[1]='j'
    elif m[1]=='k0':idx=['('+idx[0]+')+i1',idx[1]]
    else:idx=['('+idx[0]+')+i1']+['('+x+')+1' for x in idx[1:]]
   else:idx=['native_i' if x==':' else x for x in idx]
   ref=(m[1],tuple(poly(x) for x in idx));token='array_'+str(len(refs));refs.append(ref);s=s[:m.start()]+token+s[end:]
  if key:
   assert len(refs)==1,(name,s,refs)
   return refs[0]
  node=ast.parse(s,mode='eval').body
  def norm(n):
   if isinstance(n,ast.Name):return refs[int(n.id[6:])] if n.id.startswith('array_') else ('name',n.id)
   if isinstance(n,ast.Constant):return ('constant',float(n.value).hex())
   if isinstance(n,ast.UnaryOp):return (type(n.op).__name__,norm(n.operand))
   if isinstance(n,ast.BinOp):
    if isinstance(n.op,ast.Pow) and isinstance(n.right,ast.Constant) and n.right.value==2:return ('Mult',norm(n.left),norm(n.left))
    return (type(n.op).__name__,norm(n.left),norm(n.right))
   if isinstance(n,ast.Call):
    if n.func.id=='real_t' and len(n.args)==1 and isinstance(n.args[0],ast.Name) and n.args[0].id in {'zero','one_half'}:return norm(n.args[0])
    args=[]
    for a in n.args:
     v=norm(a)
     if n.func.id in {'min','max'} and v[:2]==('call',n.func.id):args+=list(v[2])
     else:args.append(v)
    return ('call',n.func.id,tuple(args))
   raise ValueError(ast.dump(n))
  return norm(node)
 native=[]
 for st in statements(NB[name]):
  m=re.match(r'^(\w+)(\(.*?\))?\s*=(?!=)\s*(.*)$',st)
  if not m:continue
  var,ind,rhs=m.groups()
  if var in {'i','k','kk','kl','km1','it','lmt','k0found'} or var=='k0':continue
  if var not in meta:continue
  lhs=expr(var+ind,key=True) if ind else ('scalar',var)
  native.append((lhs,expr(rhs),st))
 observed=[]
 # Capture initialized scalar declarations, plain assignments and array stores.
 names=set(meta)-{'i','k','kk','kl','km1','it','lmt','k0found','k0'}
 for m in re.finditer(r'(?<![\w.])('+'|'.join(sorted(names,key=len,reverse=True))+r')\b',CB[name]):
  var=m[1];p=m.end()
  while CB[name][p].isspace():p+=1
  index=''
  if CB[name][p]=='(':
   raw,p=balanced(CB[name],p);index='('+raw+')'
  while CB[name][p].isspace():p+=1
  if CB[name][p]!='=' or CB[name][p:p+2]=='==':continue
  rhs=CB[name][p+1:CB[name].index(';',p)].strip()
  # Adjacent simple declarations d1=...,d2=... retain two science assignments.
  rhs=split_top(rhs)[0]
  # map1_ppm's input kernel uses zero-based k; later science kernels retain Native k.
  if name=='map1_ppm' and m.start()<CB[name].index('ppm2m('):
   index=re.sub(r'\bk\b','(k-1)',index);rhs=re.sub(r'\bk\b','(k-1)',rhs)
  try:lhs=expr(var+index,True,True) if index else ('scalar',var);tree=expr(rhs,True)
  except (SyntaxError,ValueError):continue
  observed.append((lhs,tree,rhs))
 expected={};actual={}
 for lhs,rhs,raw in native:expected.setdefault(lhs,[]).append((rhs,raw))
 for lhs,rhs,raw in observed:actual.setdefault(lhs,[]).append((rhs,raw))
 mismatches=[]
 for key,rows in expected.items():
  got=actual.get(key,[])
  if [r[0] for r in rows]!=[r[0] for r in got]:mismatches.append({'lhs':repr(key),'native':[r[1] for r in rows],'cpp':[r[1] for r in got]})
 reports[name]={'Native_scientific_assignment_trees':len(native),'Cpp_scientific_assignment_trees':sum(len(actual.get(k,[])) for k in expected),'mismatches':mismatches}
record={'provenance':'Owner-written static expression/index check; not independent numerical evidence.','cpp_sha256':hashlib.sha256((D/'remapping_module.cpp').read_bytes()).hexdigest(),'routines':reports}
(D/'formula_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2));assert not any(p['mismatches'] for p in reports.values())
