"""Prepare the isolated Native-order C++ draft; source audit is a separate step."""
from pathlib import Path
from dataclasses import dataclass,field
import re,sys,json,textwrap
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import declaration,statements,find_procedure
from generate_derived_types import split_top
source=(root/'src/CLUBB_core/diffusion.F90').read_text()
names=re.findall(r'^  subroutine (\w+)',source,re.M)
blocks={n:re.search(r'^  subroutine '+n+r'\s*\(.*?^  end subroutine '+n+r'\b',source,re.M|re.S).group() for n in names}

@dataclass
class Node:
 kind:str
 text:str=''
 body:list=field(default_factory=list)
 other:list=field(default_factory=list)

# Keep comments and source continuation order, including conditional writes.
def tokens(text):
 out=[];pending='';comments=[];in_body=False
 for line in text.splitlines():
  stripped=line.strip()
  if 'begin code' in stripped.lower():in_body=True
  if stripped.startswith('!$acc'):continue
  if stripped.startswith('!'):
   if pending:comments.append(Node('comment',stripped[1:]));continue
   out.append(Node('comment',stripped[1:]));continue
  if not stripped:
   if not pending:out.append(Node('blank'))
   continue
  quote='';cut=len(line)
  for j,c in enumerate(line):
   if c in "'\"":quote='' if quote==c else c if not quote else quote
   if c=='!' and not quote:cut=j;break
  code=line[:cut].strip()
  if not code:continue
  pending+=' '+code.strip('& ')
  if in_body and cut<len(line):
   note=line[cut+1:].strip()
   if note and not re.match(r'(?:intent\s*\(|(?:Optional\s+)?In(?:Out)?\b|Out\b)',note,re.I):comments.append(Node('comment',' '+note))
  if code.endswith('&'):continue
  if re.match(r'(?:else|case)\b',pending.strip(),re.I):
   out.append(Node('stmt',pending.strip()));out.extend(comments)
  else:
   out.extend(comments);out.append(Node('stmt',pending.strip()))
  comments=[];pending=''
 assert not pending
 return out

def parse(ts,pos=0,stops=(),single=False):
 out=[]
 while pos<len(ts):
  t=ts[pos];s=t.text;low=s.lower()
  if t.kind in {'comment','blank'}:out.append(t);pos+=1;continue
  if any(re.match(stop,low) for stop in stops):return out,pos
  if re.match(r'if\s*\(',s,re.I):
   start=s.index('(');d=1;i=start+1
   while d:d+=(s[i]=='(')-(s[i]==')');i+=1
   cond=s[start+1:i-1];tail=s[i:].strip()
   if tail.lower()=='then':
    body,pos=parse(ts,pos+1,(r'else\b',r'else\s*if',r'elseif',r'end\s*if'))
    other=[]
    if re.match(r'else\s*if|elseif',ts[pos].text,re.I):
     ts[pos].text=re.sub(r'^else\s*if|^elseif','if',ts[pos].text,flags=re.I)
     nested,pos=parse(ts,pos,single=True)
     other=nested;pos-=1
    elif ts[pos].text.lower()=='else':other,pos=parse(ts,pos+1,(r'end\s*if',))
    assert re.match(r'end\s*if',ts[pos].text,re.I),(s,ts[pos].text)
    out.append(Node('if',cond,body,other));pos+=1
    if single:return out,pos
    continue
   out.append(Node('if',cond,[Node('stmt',tail)]));pos+=1;continue
  if re.match(r'do\s+\w+\s*=',s,re.I):
   body,pos=parse(ts,pos+1,(r'end\s*do',));out.append(Node('for',s,body));pos+=1;continue
  if re.match(r'select\s+case',s,re.I):
   cases=[];pos+=1
   while not re.match(r'end\s*select',ts[pos].text,re.I):
    if ts[pos].kind in {'comment','blank'}:
     if ts[pos].kind=='comment':cases.append(ts[pos])
     pos+=1;continue
    assert re.match(r'case',ts[pos].text,re.I),ts[pos]
    case=ts[pos].text;body,pos=parse(ts,pos+1,(r'case',r'end\s*select'))
    cases.append(Node('case',case,body))
   out.append(Node('select',s,cases));pos+=1;continue
  out.append(t);pos+=1
 return out,pos

metadata={};args={};bodies={};headers={};initializers={};variable_docs={}
for n,b in blocks.items():
 ts=tokens(b);signature=next(t.text for t in ts if t.kind=='stmt');args[n]=split_top(signature[signature.index('(')+1:signature.rindex(')')])
 lastdecl=max(i for i,t in enumerate(ts) if t.kind=='stmt' and '::' in t.text)
 info={};init={}
 for t in ts[:lastdecl+1]:
  if t.kind=='stmt' and '::' in t.text:
   for key,a in declaration(t.text).items():
    info[key]=a
   if re.search(r'\bparameter\b',t.text,re.I):
    for item in split_top(t.text.split('::')[1]):
     key,val=item.split('=',1);init[key.strip().lower()]=val.strip()
 # This source declares INTENT for every argument; no contract overrides.
 metadata[n]=info;initializers[n]=init
 variable_docs[n]={}
 last_marker=b.rfind('::');blank=b.find('\n\n',last_marker)
 documentation_prefix=b[b.lower().index('implicit none'):blank if blank>=0 else len(b)]
 for key,a in info.items():
  for raw in documentation_prefix.splitlines():
   if '!' not in raw or raw.strip().startswith('!'):continue
   code,comment=raw.split('!',1)
   if re.search(r'(?:\b|::\s*)'+re.escape(a.name)+r'\s*(?:,\s*&|&)?\s*$',code):
    variable_docs[n][key]=comment.strip();break
 headers[n]=[t for t in ts[1:lastdecl+1] if t.kind=='comment']
 bodies[n],end=parse(ts[lastdecl+1:-1]);assert end==len(ts[lastdecl+1:-1]),n

# Keep declaration-comment continuations attached to their documented argument.
for n,block in blocks.items():
 lines=block[block.lower().index('implicit none'):].splitlines();remove=[]
 for key,a in metadata[n].items():
  if key not in variable_docs[n]:continue
  for i,line in enumerate(lines):
   if '!' not in line or line.strip().startswith('!'):continue
   code,comment=line.split('!',1)
   if not re.search(r'(?:\b|::\s*)'+re.escape(a.name)+r'\s*(?:,\s*&|&)?\s*$',code):continue
   following=[]
   for raw in lines[i+1:]:
    if not raw.strip().startswith('!') or raw.strip().startswith('!$acc'):break
    content=raw.strip()[1:]
    if not content.strip() or re.match(r'[-=]{3}',content.strip()):break
    following.append(content)
   if following:
    variable_docs[n][key]+=' '+' '.join(part.strip() for part in following)
    remove.extend(following)
   break
 for content in remove:
  hit=next((j for j,t in enumerate(headers[n]) if t.text==content),None)
  if hit is not None:headers[n].pop(hit)

# Retain Native local parameter labels beside their exact definitions.
for n,block in blocks.items():
 for key in initializers[n]:
  for line in block.splitlines():
   if '!' not in line or line.strip().startswith('!'):continue
   code,comment=line.split('!',1)
   if re.match(r'\s*'+re.escape(key)+r'\s*=',code,re.I):
    variable_docs[n][key]=comment.strip();break

variable_docs['diffusion_zm_lhs'].update({'k_zm':'Coef. of eddy diffusivity at momentum levels [m^2/s]','invrs_rho_ds_zm':'Inverse dry, static density on momentum levels [m^3/kg]','rho_ds_zt':'Dry, static density on thermodynamic levels [kg/m^3]'})
# All original named state owners stay external; metadata is used only to map
# their component arrays without capturing an owning UDT in a device lambda.
component_ranks={'err_info.err_code':1,'err_info.err_header':2,'gr.zm':2,'gr.zt':2,
 'gr.invrs_dzm':2,'gr.invrs_dzt':2,'gr.dzm':2,'gr.dzt':2,'gr.weights_zm2zt':3,'gr.weights_zt2zm':3}
for fieldname in re.findall(r'nu_vert_res_dep%(\w+)',source):component_ranks['nu_vert_res_dep.'+fieldname]=1
module_constants={}
parameter_indices=set(re.findall(r'inline constexpr int (i\w+)',(root/'clubb_kokkos/src/CLUBB_core/parameter_indices.hpp').read_text()))

# Power lowering preserves explicit source grouping; no FMA/reassociation.
def powers(s):
 while '**' in s:
  p=s.rindex('**');end=p
  while s[end-1].isspace():end-=1
  start=end-1
  if s[start]==')':
   dep=1;start-=1
   while dep:dep+=(s[start]==')')-(s[start]=='(');start-=1
   start+=1;j=start-1
   while j>=0 and (s[j].isalnum() or s[j] in '_.%'):j-=1
   if j<start-1:start=j+1
  else:
   while start>=0 and (s[start].isalnum() or s[start] in '_.%'):start-=1
   start+=1
  e=p+2
  while s[e].isspace():e+=1
  begin=e
  if s[e]=='(':
   dep=1;e+=1
   while dep:dep+=(s[e]=='(')-(s[e]==')');e+=1
  else:
   while e<len(s) and (s[e].isalnum() or s[e] in '_.'):e+=1
  base=s[start:end];ex=s[begin:e]
  replacement='('+base+' * '+base+')' if ex=='2' else '(('+base+' * '+base+') * '+base+')' if ex=='3' else 'Kokkos::pow('+base+','+ex+')'
  s=s[:start]+replacement+s[e:]
 return s

def split_call(s,start):
 d=1;i=start
 while d:d+=(s[i]=='(')-(s[i]==')');i+=1
 return split_top(s[start:i-1]),i

def split_payload(s):
 quote='';depth=0;start=0;out=[]
 for i,ch in enumerate(s):
  if ch in "\"'":quote='' if quote==ch else ch if not quote else quote
  if quote:continue
  if ch=='(':depth+=1
  elif ch==')':depth-=1
  elif ch==',' and depth==0:out.append(s[start:i].strip());start=i+1
 out.append(s[start:].strip());return out

def format_call(name,values,intents):
 labels={'in':'In','inout':'InOut','out':'Out'}
 groups=[];current=[];last=None
 for value,intent in zip(values,intents):
  if current and (intent!=last or len(', '.join(current+[value]))>80):groups.append((current,last));current=[]
  current.append(value);last=intent
 if current:groups.append((current,last))
 return name+'(\n'+'\n'.join(', '.join(v)+(',' if i<len(groups)-1 else ');')+' // '+labels[t] for i,(v,t) in enumerate(groups))

class Emitter:
 def __init__(self,n):
  self.name=n;self.meta=metadata[n];self.counter=0;self.in_kernel=False;self.coords=None;self.aliases={};self.host_aliases={}
  executable='\n'.join(statements(blocks[n]))
  for obj,fieldname in re.findall(r'(gr|pdf_implicit_coefs_terms|pdf_params|sclr_idx|err_info|nu_vert_res_dep)%(\w+)',executable):self.aliases[obj+'.'+fieldname]=obj+'_'+fieldname
 def rank(self,var):
  return len(self.meta[var.lower()].dims) if var.lower() in self.meta else component_ranks.get(var,0)
 def expr(self,s,host=False,coords=None):
  s=s.strip().replace('&',' ').replace('%','.')
  s=re.sub(r'\bcp\b','Cp',s)
  # String literals are carried without rewriting identifiers within them.
  strings=[]
  def stash(m):strings.append(json.dumps(m[0][1:-1].replace("''","'")));return '@STR'+str(len(strings)-1)+'@'
  s=re.sub(r'"[^"\n]*"|\'[^\'\n]*\'',stash,s)
  s=powers(s)
  s=re.sub(r'((?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)_core_rknd',lambda m:'real_t('+m[1]+')',s)
  for old,new in [('.and.','&&'),('.or.','||'),('.not.','!'),('.true.','true'),('.false.','false'),('/=','!=')]:s=re.sub(re.escape(old),new,s,flags=re.I)
  s=re.sub(r'\bsign\s*\(', 'Kokkos::copysign(', s, flags=re.I)
  for fun in ['exp','abs','sqrt','log','log10']:s=re.sub(r'(?<![\w:])'+fun+r'\s*\(','Kokkos::'+fun+'(',s,flags=re.I)
  # Intrinsic lowering before View indexing handles nested function operands.
  for fun,op in [('min','<'),('max','>')]:
   while m:=re.search(r'(?<![\w:])'+fun+r'\s*\(',s,re.I):
    items,end=split_call(s,m.end());assert len(items)==2,s
    a,b=items;s=s[:m.start()]+'(('+a+') '+op+' ('+b+') ? ('+a+') : ('+b+'))'+s[end:]
  s=re.sub(r'len_trim\(\s*(\w+)\s*\)',r'\1.size()',s,flags=re.I)
  s=re.sub(r'present\s*\(\s*(\w+)\s*\)',r'(\1 != nullptr)',s,flags=re.I)
  s=re.sub(r'epsilon\(\s*rtm\([^()]*\)\s*\)', 'std::numeric_limits<real_t>::epsilon()',s,flags=re.I)
  s=s.replace('stats.l_sample','stats.flags().l_sample')
  s=re.sub(r'clubb_at_least_debug_level_api\s*\(',r'f2k::clubb_at_least_debug_level_api(',s)
  s=re.sub(r'\bvar_on_stats_list\s*\(',r'f2k::var_on_stats_list(',s)
  # Native ANY is a host diagnostic/reduction; existing status owner handles equality.
  s=re.sub(r'any\(\s*err_info.err_code\s*==\s*clubb_fatal_error\s*\)',r'any_equal(err_info.err_code, clubb_fatal_error)',s,flags=re.I)
  for pattern in [r'any\(\s*(xm_old|xm)\s*(?:\(:,\s*:\))?\s*<\s*(zero|xm_threshold)\s*\)']:
   s=re.sub(pattern,lambda m:'any_less('+m[1]+','+m[2]+')',s,flags=re.I)
  # Arrays, including component Views; ordinary function calls remain intact.
  allvars={a.name:len(a.dims) for a in self.meta.values() if a.dims}|component_ranks
  for var,rank in sorted(allvars.items(),key=lambda x:-len(x[0])):
   pattern=r'(?<![\w.])'+re.escape(var)+r'\s*\('
   while m:=re.search(pattern,s):
    indices,end=split_call(s,m.end());assert len(indices)==rank,(self.name,var,indices,rank)
    base=self.aliases.get(var,var) if self.in_kernel else var
    if any(':' in index for index in indices):
     if coords is not None:
      generated=[];j=0
      for axis,idx in enumerate(indices):
       if ':' in idx:
        parts=idx.split(':');lo=parts[0].strip() or '1';step=parts[2].strip() if len(parts)==3 else '1'
        generated.append('('+self.expr(lo)+' - 1 + '+coords[j]+'*('+self.expr(step)+'))');j+=1
       else:generated.append('('+self.expr(idx)+' - 1)' if not idx.strip() in parameter_indices else self.expr(idx))
      value=base+'('+', '.join(generated)+')'
     else:
      items=[]
      for axis,idx in enumerate(indices):
       if idx.strip()==':':items.append('Kokkos::ALL')
       elif ':' in idx:
        parts=idx.split(':');assert len(parts)==2,('stride subview',self.name,var,idx)
        items.append('std::make_pair('+self.expr(parts[0] or '1')+'-1,'+self.expr(parts[1] or str(self.meta[var.lower()].dims[axis]))+')')
       else:items.append('('+self.expr(idx)+'-1)' if not idx.strip() in parameter_indices else self.expr(idx))
      value='Kokkos::subview('+base+','+','.join(items)+')'
    else:
     value=base+'('+', '.join(self.expr(idx)+( '' if idx.strip() in parameter_indices else ' - 1') for idx in indices)+')'
     if host and rank:value=self.host_aliases.get(base,'to_host_space('+base+')')+value[len(base):]
    # Mark this converted reference so the same name is not repeatedly matched.
    marker='@ARR'+str(len(strings))+'@';strings.append(value);s=s[:m.start()]+marker+s[end:]
   if coords is not None:
    replacement=(self.aliases.get(var,var) if self.in_kernel else var)+'('+', '.join(coords[:rank])+')'
    s=re.sub(r'(?<![\w.])'+re.escape(var)+r'(?![\w.(])',lambda m:replacement,s)
  if self.in_kernel:
   for var,alias in self.aliases.items():s=re.sub(r'(?<![\w.])'+re.escape(var)+r'(?![\w.])',alias,s)
  for opt in ['wpxp_cl_num']:
   if opt in self.meta and self.meta[opt].optional:s=re.sub(r'(?<![\w])'+opt+r'(?!\w|\s*!=)',r'(*'+opt+')',s)
  s=s.replace('//','+') # Native character concatenation.
  # Restore nested placeholders in reverse order.
  for j in reversed(range(len(strings))):s=s.replace('@ARR'+str(j)+'@',strings[j]).replace('@STR'+str(j)+'@',strings[j])
  s=re.sub(r'\bthreshold_hi\b(?!\s*!=\s*nullptr)', 'threshold_hi_value',s) if 'threshold_hi' in self.meta else s
  return re.sub(r'\s+',' ',s).strip()
 def emit(self,nodes,kernel=False):
  self.in_kernel=kernel;out=[]
  for node in nodes:
   kind,s=node.kind,node.text
   if kind=='raw':out.append(s);continue
   if kind=='blank':out.append('');continue
   if kind=='comment':out.append('//'+s);continue
   if kind=='if':
    if 'sponge_damp_settings%' in s:
     out.append('// Native '+s.strip()+':\n// Sponge damping is excluded with the JAX scope.\n// TODO: Reject sponge-enabled cases upstream; see INTEGRATION.md.');continue
    cond=self.expr(s,host=not kernel)
    if not any(x.kind!='blank' for x in node.body):node.body=[Node('comment',' Native OpenACC host update: the stats API stages its own inputs.')]
    out.append('if ('+cond+') {\n'+self.emit(node.body,kernel)+'\n}')
    if node.other:
     out[-1]+= (' else '+self.emit(node.other,kernel)) if len(node.other)==1 and node.other[0].kind=='if' else ' else {\n'+self.emit(node.other,kernel)+'\n}'
    continue
   if kind=='select':
    condition=s[s.index('(')+1:s.rindex(')')];out.append('switch ('+self.expr(condition)+') {')
    for case in node.body:
     if case.kind=='comment':out.append('//'+case.text);continue
     val='default' if 'default' in case.text.lower() else 'case '+self.expr(case.text[case.text.index('(')+1:case.text.rindex(')')])
     out.append(val+': {\n'+self.emit(case.body,kernel)+'\nbreak;\n}')
    out.append('}');continue
   if kind=='for':
    m=re.match(r'do\s+(\w+)\s*=\s*(.*)',s,re.I);var=m[1];ranges=split_top(m[2]);lo,hi=ranges[:2];step=ranges[2] if len(ranges)>2 else '1'
    # Only kernels with ordinary scalar assignments/control are fused here.
    def ordinary(ns):
     for x in ns:
      if x.kind=='stmt' and (re.match(r'(?:call|write|return)\b',x.text,re.I) or ':' in x.text):return False
      if x.kind in {'for','if','select','case'} and not ordinary(x.body+x.other):return False
     return True
    if not kernel and ordinary(node.body):
     loops=[(var,lo,hi,step)];body=node.body
     while len([x for x in body if x.kind not in {'comment','blank'}])==1 and next(x for x in body if x.kind not in {'comment','blank'}).kind=='for':
      nested=next(x for x in body if x.kind not in {'comment','blank'});m=re.match(r'do\s+(\w+)\s*=\s*(.*)',nested.text,re.I);rs=split_top(m[2]);loops.append((m[1],rs[0],rs[1],rs[2] if len(rs)>2 else '1'));body=nested.body
     self.counter+=1;idx='cell';self.in_kernel=True
     bounds=[]
     for _,a,b,c in loops:
      if a=='1' and c=='1':bounds.append(self.expr(b))
      elif a=='2' and b=='nzm-1' and c=='1':bounds.append('nzm-2')
      else:bounds.append(f'(({self.expr(b)} - ({self.expr(a)})) / ({self.expr(c)}) + 1)')
     count=' * '.join('std::max(0,'+bound+')' for bound in bounds);parts=['Kokkos::parallel_for("'+self.name+'::loop'+str(self.counter)+'", '+count+', KOKKOS_LAMBDA(const int cell) {','int remaining=cell;']
     for (v,a,b,c),extent in reversed(list(zip(loops,bounds))):parts.append('const int '+v+'='+self.expr(a)+'+remaining%('+extent+')'+('' if c=='1' else '*('+self.expr(c)+')')+'; remaining/=('+extent+');')
     # Private scalar temporaries in native loops cannot mutate const captures.
     text='\n'.join(x.text for x in walk(body));loopvars={v for v,*_ in loops}
     for a in self.meta.values():
      if not a.dims and a.name not in args[self.name] and a.name not in loopvars and re.search(r'\b'+a.name+r'\s*=',text) and a.name.lower() not in initializers[self.name]:parts.append(scalar_type(a)+' '+a.name+';')
     parts.append(self.emit(body,True));parts.append('});');out.append('\n'.join(parts));self.in_kernel=kernel;continue
    cond=var+(' <= ' if not step.startswith('-') else ' >= ')+self.expr(hi) if re.fullmatch(r'-?\d+',step) else '('+self.expr(step)+'>0 ? '+var+' <= '+self.expr(hi)+' : '+var+' >= '+self.expr(hi)+')'
    out.append('for ('+var+'='+self.expr(lo)+'; '+cond+'; '+var+' += '+self.expr(step)+') {\n'+self.emit(node.body,kernel)+'\n}');continue
   if kind=='stmt':out.append(self.statement(s,kernel));self.in_kernel=kernel;continue
   raise ValueError((kind,s))
  return '\n'.join(out)
 def statement(self,s,kernel):
  if s.lower()=='return':return 'return;'
  if re.match(r'write\s*\(',s,re.I):
   start=s.index('(');items,end=split_call(s,start+1);payload=split_payload(s[end:]);parts=[]
   for item in payload:
    if re.match(r'new_line\(',item,re.I):parts.append('"\\n"');continue
    if re.match(r'err_info%err_header\(',item,re.I):
     index=item[item.index('(')+1:item.rindex(')')]
     parts.append('dump_header(err_info.err_header,'+self.expr(index)+'-1)');continue
    converted=self.expr(item,host=True)
    # Host dumps preserve native View element ordering through one print adapter.
    if re.match(r'^[\w.]+$',converted) and self.rank(converted):parts.append('dump_values('+converted+')')
    elif converted.startswith('Kokkos::subview('):parts.append('dump_values('+converted+')')
    elif converted.startswith('err_info.err_header('):parts.append('dump_header(err_info.err_header,'+str(item)+')')
    else:parts.append(converted)
   return '{ std::ostringstream diagnostic_message; diagnostic_message<<std::setprecision(std::numeric_limits<real_t>::max_digits10);\ndiagnostic_message << '+ ' << " " << '.join(parts)+' << "\\n";\nstd::cerr<<diagnostic_message.str(); }'
  if re.match(r'call\s+',s,re.I):
   name=s.split()[1].split('(')[0];a,end=split_call(s,s.index('(')+1)
   if name=='band_solve':raise ValueError('special matrix boundary')
   converted=[]
   if name in names:
    contract=[metadata[name][key.lower()] for key in args[name]]
    positionals={};p=0
    for val in a:
     if re.match(r'^\w+\s*=(?!=)',val):key,value=val.split('=',1);positionals[key.strip().lower()]=value
     else:positionals[contract[p].name.lower()]=val;p+=1
    for param in contract:
     value=positionals.get(param.name.lower())
     if value is None:assert param.optional;converted.append('nullptr');continue
     expr=self.expr(value,host=True)
     if param.optional:
      expr=expr if value.strip() in self.meta and self.meta[value.strip()].optional else '&'+expr
     converted.append(expr)
   else:
    converted=[self.expr(val,host=True) for val in a]
   if name.startswith('stats_'):
    # Array arithmetic is evaluated visibly before crossing the stats bridge.
    raw=a[1]
    if re.search(r'[/+*-]',raw):
     matches=[key for key,val in self.meta.items() if val.dims and re.search(r'\b'+re.escape(val.name)+r'\b',raw)]
     array=self.meta[matches[0]].name;rank=self.rank(array);assert rank in {1,2},(s,array,rank)
     self.counter+=1;label='stats_values_'+str(self.counter);self.in_kernel=True
     indices=['i','k'][:rank];value=self.expr(raw,coords=indices);self.in_kernel=False
     dims=','.join(array+'.extent('+str(axis)+')' for axis in range(rank))
     indexcode='const int i=cell;' if rank==1 else 'const int i=cell%'+array+'.extent(0); const int k=cell/'+array+'.extent(0);'
     values='f2k::Real1D('+label+')' if rank==1 else label
     return '{\nReal'+str(rank)+'D '+label+'(Kokkos::view_alloc(Kokkos::WithoutInitializing,"'+self.name+'::'+label+'"), '+dims+');\nKokkos::parallel_for("'+self.name+'::stats_value", '+array+'.size(), KOKKOS_LAMBDA(int cell) { '+indexcode+' '+label+'('+','.join(indices)+')='+value+'; });\nf2k::'+name+'('+converted[0]+','+values+', // In\n'+converted[2]+'); // InOut\n}'
   if name in names:intents=[metadata[name][key.lower()].intent for key in args[name]]
   elif name.startswith('stats_'):intents=['in','in','inout']
   else:
    owner={'calc_xpwp':'calc_xpwp_2D'}.get(name,name);intents=[a.intent for a in find_procedure(owner).arguments]
   if name in public_lower or (name.startswith('stats_') and name not in names):name='f2k::'+name
   return format_call(name,converted,intents)
  m=re.match(r'([\w%.]+)(\(.*?\))?\s*=\s*(.*)',s,re.S)
  if m:
   var=m[1].replace('%','.');index=m[2];value=m[3];rank=self.rank(var)
   # Array-valued whole assignment or selected native array sections.
   if rank and (index is None or ':' in index):
    left=self.expr(var+(index or ''))
    if re.match(r'\s*(zt2zm_api|zm2zt_api|ddzt|ddzm|zt2zm2zt|zm2zt2zm)\s*\(',value):
     value=self.expr(value)
     for fun in []:
      if value.startswith(fun+'('):
       aa,_=split_call(value,value.index('(')+1);aa[-1]='to_execution_space('+aa[-1]+')';value=fun+'('+','.join(aa)+')'
     return 'Kokkos::deep_copy('+left+','+value+');'
    if var=='H_zagl' and value.startswith('smooth_heaviside_peskin'):
     return 'Kokkos::deep_copy('+left+','+self.expr(value)+');'
    if var=='err_info.err_code':return 'Kokkos::deep_copy(err_info.err_code,'+self.expr(value)+');'
    # Preserve aliasing semantics of Native reversed whole-array assignment.
    if index and ':-1' in value:raise ValueError('special reverse assignment')
    selected=split_top(index[1:-1]) if index else [':']*rank
    dimensions=[]
    for axis,item in enumerate(selected):
     if ':' in item:
      lo,*rest=item.split(':');hi=rest[0] if rest and rest[0].strip() else var+'.extent('+str(axis)+')';lo=lo or '1';step=rest[1] if len(rest)>1 else '1'
      dimensions.append('('+self.expr(hi)+'-('+self.expr(lo)+'))/('+self.expr(step)+')+1')
    coords=['slice_'+str(j) for j in range(len(dimensions))]
    self.counter+=1;self.in_kernel=True
    lv=self.expr(var+(index or ''),coords=coords);rv=self.expr(value,coords=coords)
    self.in_kernel=False
    parts=['Kokkos::parallel_for("'+self.name+'::array'+str(self.counter)+'", '+' * '.join('('+x+')' for x in dimensions)+', KOKKOS_LAMBDA(int cell) {','int remaining=cell;']
    for c,d in zip(coords,dimensions):parts.append('const int '+c+'=remaining%('+d+'); remaining/=('+d+');')
    parts.append(lv+' = '+rv+';');parts.append('});');return '\n'.join(parts)
   if rank and not kernel:self.in_kernel=True
   left=self.expr(var+(index or ''));right=self.expr(value)
   self.in_kernel=kernel
   if rank and not kernel:
    return 'Kokkos::parallel_for("'+self.name+'::set",1,KOKKOS_LAMBDA(int) { '+left+'='+right+'; });'
   return left+' = '+right+';'
  raise ValueError((self.name,s))

def walk(nodes):
 for n in nodes:
  yield n
  yield from walk(n.body);yield from walk(n.other)
def scalar_type(a):return {'real':'real_t','double':'double','integer':'int','logical':'bool','character':'std::string'}.get(a.kind,a.kind)
def signature(n,prototype=False):
 arrays=[key for key in args[n] if metadata[n][key.lower()].dims and not metadata[n][key.lower()].optional];typed=[]
 for key in args[n]:
  a=metadata[n][key.lower()]
  if a.dims:
   typ=('Real'+str(len(a.dims))+'D') if True else 'Array'+str(arrays.index(key))
   if a.kind=='integer':typ='Kokkos::View<int'+('*'*len(a.dims))+',Kokkos::LayoutLeft>' if True else typ
   if a.optional:typ+='*'
  elif a.kind in {'grid','sclr_idx_type','pdf_parameter','implicit_coefs_terms','nu_vertical_res_dep','err_info_type','stats_type'}:
   typ=('f2k::' if a.kind=='stats_type' else '')+a.kind+'&'
   if a.intent=='in':typ='const '+typ
  else:
   typ=('const std::string&' if a.kind=='character' else ('const ' if a.optional and a.intent=='in' else '')+scalar_type(a)+('*' if a.optional else '&' if a.intent in {'inout','out'} else ''))
  typed.append(typ+' '+a.name+(' = nullptr' if a.optional and prototype else ''))
 prefix=''

 # Preserve Native continuation groups as far as C++ typed declarations fit.
 groups=[]
 raw_signature=blocks[n].split('\n\n',1)[0]
 for line in raw_signature.splitlines():
  code=line.split('!',1)[0].strip()
  code=re.sub(r'^subroutine\s+\w+\s*\(','',code,flags=re.I)
  code=code.strip('& ').rstrip(')').strip()
  members=[x.strip() for x in code.split(',') if x.strip()]
  if members:groups.append(members)
  if ')' in line.split('!',1)[0]:break
 assert [key for group in groups for key in group]==args[n],(n,groups)
 declarations=dict(zip(args[n],typed));formatted=[];completed=0
 for group in groups:
  pending=[];intent=None
  for key in group:
   current=('optional ' if metadata[n][key.lower()].optional else '')+metadata[n][key.lower()].intent
   if pending and (current!=intent or len(', '.join(pending+[declarations[key]]))>78):
    formatted.append(', '.join(pending)+', // '+{'in':'In','inout':'InOut','out':'Out','optional in':'Optional In'}[intent]);pending=[]
   pending.append(declarations[key]);intent=current;completed+=1
  formatted.append(', '.join(pending)+(',' if completed<len(args[n]) else ')')+' // '+{'in':'In','inout':'InOut','out':'Out','optional in':'Optional In'}[intent])
 descriptions=[]
 for key in args[n]:
  note=variable_docs[n].get(key.lower(),'')
  if note and not re.fullmatch(r'(?:intent\s*\(\s*)?(?:inout|in|out)\s*\)?',note,re.I):
   descriptions.extend('// '+line for line in textwrap.wrap(key+': '+note,width=94))
 return '\n'.join(descriptions)+('\n' if descriptions else '')+prefix+'void '+n+'(\n'+'\n'.join(formatted)+'\n'


# The source operator boundary retains Native bands first. No matrix transpose.
public_lower=set()
opening=source[source.index('  ! Description:'):source.index('  implicit none')]
support='''// Compute Native left-hand-side operator coefficients.
//
// Physical inputs and scratch use column-first LayoutLeft Views. The complete
// Out operator keeps Native (band,column,level) order, with three bands ordered
// super/main/sub. Keep any solver-layout transpose at the later solver boundary.
'''+ '\n'.join('//'+line.strip()[1:] for line in opening.splitlines() if line.strip().startswith('!'))+'\n'+f'''#include "{folder.name}.hpp"
#include "constants_clubb.hpp"
'''+('#include "model_flags.hpp" // Selected l_upwind_Kh_dp_term owner remains pending.\n' if folder.name=='diffusion' else '')+'''#include <algorithm>
namespace clubb_kokkos {
'''
if module_constants:
 support+='namespace {\n'
 for key,value in module_constants.items():support+='constexpr int '+key+'='+str(value)+';\n'
 support+='} // namespace\n'
for n in names:
 e=Emitter(n);support+='//=============================================================================\n'+signature(n)+' {\n'
 support+='\n'.join('//'+t.text for t in headers[n])+'\n'
 for key,a in metadata[n].items():
  if a.name in args[n]:continue
  doc=' // '+variable_docs[n][key] if key in variable_docs[n] else ''
  if a.dims:support+='Kokkos::View<'+scalar_type(a)+'*'*len(a.dims)+',Kokkos::LayoutLeft> '+a.name+'(Kokkos::view_alloc(Kokkos::WithoutInitializing,"'+n+'::'+a.name+'"),'+','.join(e.expr(dim) for dim in a.dims)+');'+doc+'\n'
  else:support+=('constexpr ' if key in initializers[n] else '')+scalar_type(a)+' '+a.name+('='+e.expr(initializers[n][key]) if key in initializers[n] else '')+';'+doc+'\n'
 body=e.emit(bodies[n])
 for key,alias in e.aliases.items():
  if re.search(r'\b'+alias+r'\b',re.sub(r'//[^\n]*','',body)):support+='const auto '+alias+'='+key+';\n'
 support+=body+'\n}\n\n'
support+='} // namespace clubb_kokkos\n'
header='''// Fill Native band operators while retaining caller-owned output Views.
//
// diffusion_zt_lhs and diffusion_zm_lhs fill caller-owned lhs arrays in Native
// (band,column,level) order, with super/main/sub bands. Physical diffusivity,
// density and spacing fields are column first. All arguments are required. Generic
// entry overloads preserve existing caller sections through the shared array
// owner, without any band transpose. This draft has no build/runtime evidence;
// see INTEGRATION.md before promotion.
#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"
namespace clubb_kokkos {
'''
for n in names:header+=signature(n,True)+';\n\n'
header+='// Preserve actual operator-call sections at the public array API boundary.\n'
for n in names:
 arrays=[metadata[n][key.lower()] for key in args[n] if metadata[n][key.lower()].dims]
 types=[]
 for key in args[n]:
  a=metadata[n][key.lower()]
  typ='Array'+str(arrays.index(a)) if a.dims else 'const grid&' if a.kind=='grid' else scalar_type(a)
  types.append(typ+' '+a.name)
 header+='template<'+','.join('class Array'+str(j) for j in range(len(arrays)))+'>\nvoid '+n+'('+','.join(types)+') {\n'
 for a in arrays:header+='auto '+a.name+'_exec=to_execution_space('+a.name+(',false' if a.intent=='out' else '')+');\n'
 header+=n+'('+','.join(key+'_exec' if metadata[n][key.lower()].dims else key for key in args[n])+');\n'
 for a in arrays:
  if a.intent!='in':header+='from_execution_space('+a.name+','+a.name+'_exec);\n'
 header+='}\n'
header+='} // namespace clubb_kokkos\n'
(folder/(folder.name+'.cpp')).write_text(support)
(folder/(folder.name+'.hpp')).write_text(header)
print('Prepared both complete Native '+folder.name+' operators in isolated draft.')
