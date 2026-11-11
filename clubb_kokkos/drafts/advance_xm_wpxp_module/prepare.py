"""Prepare the isolated Native-order C++ draft; source audit is a separate step."""
from pathlib import Path
from dataclasses import dataclass,field
import re,sys,json,textwrap
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import declaration,statements,find_procedure
from generate_derived_types import split_top
source=(root/'src/CLUBB_core/advance_xm_wpxp_module.F90').read_text()
names=re.findall(r'^  subroutine (\w+)',source,re.M)
blocks={n:re.search(r'^  subroutine '+n+r'\(.*?^  end subroutine '+n+r'\b',source,re.M|re.S).group() for n in names}

@dataclass
class Node:
 kind:str
 text:str=''
 body:list=field(default_factory=list)
 other:list=field(default_factory=list)

# Keep comments and source continuation order, including conditional writes.
def tokens(text):
 out=[];pending='';comments=[]
 for line in text.splitlines():
  stripped=line.strip()
  if stripped.startswith('!$acc'):continue
  if stripped.startswith('!'):
   if pending:comments.append(Node('comment',stripped[1:]));continue
   out.append(Node('comment',stripped[1:]));continue
  if not stripped:continue
  quote='';cut=len(line)
  for j,c in enumerate(line):
   if c in "'\"":quote='' if quote==c else c if not quote else quote
   if c=='!' and not quote:cut=j;break
  code=line[:cut].strip()
  if not code:continue
  pending+=' '+code.strip('& ')
  if code.endswith('&'):continue
  out.extend(comments);comments=[]
  out.append(Node('stmt',pending.strip()));pending=''
 assert not pending
 return out

def parse(ts,pos=0,stops=(),single=False):
 out=[]
 while pos<len(ts):
  t=ts[pos];s=t.text;low=s.lower()
  if t.kind=='comment':out.append(t);pos+=1;continue
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
    if ts[pos].kind=='comment':cases.append(ts[pos]);pos+=1;continue
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
 if n=='damp_coefficient':info['damped_value'].intent='out'
 metadata[n]=info;initializers[n]=init
 variable_docs[n]={}
 last_marker=b.rfind('::');blank=b.find('\n\n',last_marker)
 documentation_prefix=b[:blank if blank>=0 else len(b)]
 for key,a in info.items():
  for raw in documentation_prefix.splitlines():
   if '!' not in raw or raw.strip().startswith('!'):continue
   code,comment=raw.split('!',1)
   if re.search(r'(?:\b|::\s*)'+re.escape(a.name)+r'\s*(?:,\s*&|&)?\s*$',code):
    variable_docs[n][key]=comment.strip();break
 headers[n]=[t for t in ts[1:lastdecl+1] if t.kind=='comment']
 bodies[n],end=parse(ts[lastdecl+1:-1]);assert end==len(ts[lastdecl+1:-1]),n

# All original named state owners stay external; metadata is used only to map
# their component arrays without capturing an owning UDT in a device lambda.
component_ranks={'err_info.err_code':1,'err_info.err_header':2,'gr.zm':2,'gr.zt':2,
 'gr.invrs_dzm':2,'gr.invrs_dzt':2,'gr.dzt':2,'gr.weights_zm2zt':3,'gr.weights_zt2zm':3}
for fieldname in re.findall(r'pdf_implicit_coefs_terms%(\w+)',source):component_ranks['pdf_implicit_coefs_terms.'+fieldname]=3 if 'sclrp' in fieldname else 2
for fieldname in re.findall(r'nu_vert_res_dep%(\w+)',source):component_ranks['nu_vert_res_dep.'+fieldname]=1
module_constants={'nsub':'2','nsup':'2','xm_wpxp_thlm':'1','xm_wpxp_rtm':'2','xm_wpxp_scalar':'3','xm_wpxp_um':'4','xm_wpxp_vm':'5','ndiags2':'2','ndiags3':'3','ndiags5':'5'}
known_native={'zt2zm_api','zm2zt_api','zt2zm2zt','zm2zt2zm','ddzt','Skx_func','band_solve'}
public_lower={'calc_wp3_on_wp2','calc_turb_adv_range','diffusion_zm_lhs','term_ma_zm_lhs','diffusion_zt_lhs','term_ma_zt_lhs','xpyp_term_ta_pdf_rhs','xpyp_term_ta_pdf_lhs','xpyp_term_ta_pdf_lhs_godunov','monotonic_turbulent_flux_limit','pos_definite_adj','fill_holes_vertical_api','clip_covar'}

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
  replacement='('+base+' * '+base+')' if ex=='2' else 'Kokkos::pow('+base+','+ex+')'
  s=s[:start]+replacement+s[e:]
 return s

def split_call(s,start):
 d=1;i=start
 while d:d+=(s[i]=='(')-(s[i]==')');i+=1
 return split_top(s[start:i-1]),i

class Emitter:
 def __init__(self,n):
  self.name=n;self.meta=metadata[n];self.counter=0;self.in_kernel=False;self.coords=None;self.aliases={};self.host_aliases={}
  executable='\n'.join(statements(blocks[n]))
  for obj,fieldname in re.findall(r'(gr|pdf_implicit_coefs_terms|err_info|nu_vert_res_dep)%(\w+)',executable):self.aliases[obj+'.'+fieldname]=obj+'_'+fieldname
 def rank(self,var):
  return len(self.meta[var.lower()].dims) if var.lower() in self.meta else component_ranks.get(var,0)
 def expr(self,s,host=False,coords=None):
  s=s.strip()
  # String literals are carried without rewriting identifiers within them.
  strings=[]
  def stash(m):strings.append(m[0]);return '@STR'+str(len(strings)-1)+'@'
  s=re.sub(r'"[^"\n]*"|\'[^\'\n]*\'',stash,s)
  s=s.replace('&',' ').replace('%','.')
  s=powers(s)
  s=re.sub(r'((?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)_core_rknd',lambda m:'real_t('+m[1]+')',s)
  for old,new in [('.and.','&&'),('.or.','||'),('.not.','!'),('.true.','true'),('.false.','false'),('/=','!=')]:s=re.sub(re.escape(old),new,s,flags=re.I)
  for fun in ['exp','abs','sqrt','log','log10']:s=re.sub(r'(?<![\w:])'+fun+r'\s*\(','Kokkos::'+fun+'(',s,flags=re.I)
  # Intrinsic lowering before View indexing handles nested function operands.
  for fun,op in [('min','<'),('max','>')]:
   while m:=re.search(r'(?<![\w:])'+fun+r'\s*\(',s,re.I):
    items,end=split_call(s,m.end());assert len(items)==2,s
    a,b=items;s=s[:m.start()]+'(('+a+') '+op+' ('+b+') ? ('+a+') : ('+b+'))'+s[end:]
  s=re.sub(r'len_trim\(\s*(\w+)\s*\)',r'\1.size()',s,flags=re.I)
  s=re.sub(r'present\(\s*(\w+)\s*\)',r'(\1 != nullptr)',s,flags=re.I)
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
       else:generated.append('('+self.expr(idx)+' - 1)' if not idx.strip().startswith(('iC','ic_','ialtitude','iwpxp')) else self.expr(idx))
      value=base+'('+', '.join(generated)+')'
     else:
      items=[]
      for axis,idx in enumerate(indices):
       if idx.strip()==':':items.append('Kokkos::ALL')
       elif ':' in idx:
        parts=idx.split(':');assert len(parts)==2,('stride subview',self.name,var,idx)
        items.append('std::make_pair('+self.expr(parts[0] or '1')+'-1,'+self.expr(parts[1] or str(self.meta[var.lower()].dims[axis]))+')')
       else:items.append('('+self.expr(idx)+'-1)' if not idx.strip().startswith(('iC','ic_','ialtitude','iwpxp')) else self.expr(idx))
      value='Kokkos::subview('+base+','+','.join(items)+')'
    else:
     value=base+'('+', '.join(self.expr(idx)+( '' if idx.strip().startswith(('iC','ic_','ialtitude','iwpxp')) else ' - 1') for idx in indices)+')'
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
  s=re.sub(r'\s+',' ',s).strip()
  # Restore nested placeholders after whitespace normalization so original
  # diagnostic labels retain their ampersands, percent signs and spacing.
  for j in reversed(range(len(strings))):s=s.replace('@ARR'+str(j)+'@',strings[j]).replace('@STR'+str(j)+'@',strings[j])
  return s
 def emit(self,nodes,kernel=False):
  self.in_kernel=kernel;out=[]
  for node in nodes:
   kind,s=node.kind,node.text
   if kind=='raw':out.append(s);continue
   if kind=='comment':out.append('//'+s);continue
   if kind=='if':
    if 'sponge_damp_settings%' in s:
     out.append('// Native '+s.strip()+':\n// Sponge damping is excluded with the JAX scope. Integration must reject\n// sponge-enabled cases upstream; see INTEGRATION.md.');continue
    cond=self.expr(s,host=not kernel);snapshot='';restore=self.host_aliases.copy()
    text='\n'.join(n.text for n in walk(node.body))
    if not kernel and self.name=='advance_xm_wpxp' and 'The C7_Skw_fnc variable is outside' in text:
     snapshot='// Snapshot once for the Native host diagnostic loop.\nauto C7_Skw_fnc_host=to_host_space(C7_Skw_fnc);\n'
     self.host_aliases['C7_Skw_fnc']='C7_Skw_fnc_host'
    if not kernel and self.name=='xm_wpxp_clipping_and_stats' and s.strip().lower().startswith('any(') and 'solve_type_str' in text:
     snapshot='// Snapshot once for the Native host diagnostic loop.\nauto xm_host=to_host_space(xm);\n'
     self.host_aliases['xm']='xm_host'
    out.append('if ('+cond+') {\n'+snapshot+self.emit(node.body,kernel)+'\n}')
    self.host_aliases=restore
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
     while len([x for x in body if x.kind!='comment'])==1 and next(x for x in body if x.kind!='comment').kind=='for':
      nested=next(x for x in body if x.kind!='comment');m=re.match(r'do\s+(\w+)\s*=\s*(.*)',nested.text,re.I);rs=split_top(m[2]);loops.append((m[1],rs[0],rs[1],rs[2] if len(rs)>2 else '1'));body=nested.body
     self.counter+=1;idx='cell'
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
    out.append('for ('+var+'='+self.expr(lo)+'; '+var+(' <= ' if not step.startswith('-') else ' >= ')+self.expr(hi)+'; '+var+' += '+self.expr(step)+') {\n'+self.emit(node.body,kernel)+'\n}');continue
   if kind=='stmt':out.append(self.statement(s,kernel));self.in_kernel=kernel;continue
   raise ValueError((kind,s))
  return '\n'.join(out)
 def statement(self,s,kernel):
  if s.lower()=='return':return 'return;'
  if re.match(r'write\s*\(',s,re.I):
   start=s.index('(');items,end=split_call(s,start+1);payload=split_top(s[end:]);parts=[]
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
   return '{ std::ostringstream message; message<<std::setprecision(std::numeric_limits<real_t>::max_digits10);\nmessage << '+ ' << " " << '.join(parts)+' << "\\n";\nstd::cerr<<message.str(); }'
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
   if name=='xpyp_term_ta_pdf_rhs' and a[-1].startswith('lhs_ta_wpsclrp('):
    return '''{ // Native sequence association: this source call writes the first
// ngrdcol*nzm elements of the rank-3 LHS scalar section, not RHS.
// Preserve that source mapping; this disabled branch leaves scalar RHS undefined.
auto lhs_scalar='''+converted[-1]+''';
Kokkos::View<real_t**,Kokkos::LayoutLeft,typename decltype(lhs_scalar)::memory_space,
 Kokkos::MemoryTraits<Kokkos::Unmanaged>> lhs_sequence(lhs_scalar.data(),ngrdcol,nzm);
f2k::xpyp_term_ta_pdf_rhs('''+', '.join(converted[:-1]+['lhs_sequence'])+''');
}'''
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
     return '{\nReal'+str(rank)+'D '+label+'(Kokkos::view_alloc(Kokkos::WithoutInitializing,"'+self.name+'::'+label+'"), '+dims+');\nKokkos::parallel_for("'+self.name+'::stats_value", '+array+'.size(), KOKKOS_LAMBDA(int cell) { '+indexcode+' '+label+'('+','.join(indices)+')='+value+'; });\nf2k::'+name+'('+converted[0]+','+values+','+converted[2]+');\n}'
   if name in public_lower or name.startswith('stats_'):name='f2k::'+name
   return name+'('+', '.join(converted)+');'
  m=re.match(r'([\w%.]+)(\(.*?\))?\s*=\s*(.*)',s,re.S)
  if m:
   var=m[1].replace('%','.');index=m[2];value=m[3];rank=self.rank(var)
   # Array-valued whole assignment or selected native array sections.
   if rank and (index is None or ':' in index):
    left=self.expr(var+(index or ''))
    if re.match(r'\s*(zt2zm_api|zm2zt_api|ddzt|zt2zm2zt|zm2zt2zm)\s*\(',value):
     value=self.expr(value)
     for fun in ['ddzt','zt2zm2zt','zm2zt2zm']:
      if value.startswith(fun+'('):
       aa,_=split_call(value,value.index('(')+1);aa[-1]='to_execution_space('+aa[-1]+')';value=fun+'('+','.join(aa)+')'
     return 'Kokkos::deep_copy('+left+','+value+');'
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
   typ=('Real'+str(len(a.dims))+'D') if n==names[0] or a.optional else 'Array'+str(arrays.index(key))
   if a.kind=='integer':typ='Kokkos::View<int'+('*'*len(a.dims))+',Kokkos::LayoutLeft>' if n==names[0] else typ
   if a.optional:typ+='*'
  elif a.kind in {'grid','implicit_coefs_terms','nu_vertical_res_dep','err_info_type','stats_type'}:
   typ=('f2k::' if a.kind=='stats_type' else '')+a.kind+'&'
   if a.intent=='in':typ='const '+typ
  else:
   typ=scalar_type(a)+('*' if a.optional else '&' if a.intent in {'inout','out'} else '')
  typed.append(typ+' '+a.name+(' = nullptr' if a.optional and prototype else ''))
 prefix='template <'+','.join('class Array'+str(j) for j in range(len(arrays)))+'>\n' if n!=names[0] else ''
 if prototype and n!=names[0]:return prefix+'void '+n+'('+','.join(typed)+')'
 formatted=[]
 for i,(key,decl) in enumerate(zip(args[n],typed)):
  intent=metadata[n][key.lower()].intent.title()
  note=variable_docs[n].get(key.lower(),'')
  comment=' // '+intent+(': '+note if note else '')
  if len(decl+comment)>85 or i==len(typed)-1:
   formatted.extend('// '+line for line in textwrap.wrap(comment[4:],width=85))
   formatted.append(decl+(',' if i<len(typed)-1 else ')'))
  else:formatted.append(decl+','+comment)
 return prefix+'void '+n+'(\n'+'\n'.join(formatted)

# Special solver keeps native bands-first storage outside its existing typed
# matrix API boundary. Numerical formulas and solver selection are unchanged.
def solver_body(e):
 return '''// Solve for xm / w'x' using the native interleaved matrix equations.
// lhs is Native (band,column,row); matrix_solver_wrapper owns
// (column,band,row), so transpose the first two axes only at this boundary.
const bool reverse=l_force_descending_solves && gr.grid_dir_indx>0;
Real3D lhs_solver(Kokkos::view_alloc(Kokkos::WithoutInitializing,"xm_wpxp_solve::lhs_solver"),ngrdcol,nsup+nsub+1,2*nzm-1);
Real3D rhs_input;
if(reverse) {
 rhs_input=Real3D(Kokkos::view_alloc(Kokkos::WithoutInitializing,"xm_wpxp_solve::rhs_input"),ngrdcol,2*nzm-1,nrhs);
 Kokkos::deep_copy(rhs_input,rhs);
}
Kokkos::parallel_for("xm_wpxp_solve::matrix_order",ngrdcol*(2*nzm-1),KOKKOS_LAMBDA(int cell) {
 const int i=cell%ngrdcol,k=cell/ngrdcol;
 const int source_k=reverse?2*nzm-2-k:k;
 for(int b=0;b<nsup+nsub+1;++b) lhs_solver(i,b,k)=lhs(reverse?nsup+nsub-b:b,i,source_k);
 if(reverse)for(int n=0;n<nrhs;++n) rhs(i,k,n)=rhs_input(i,source_k,n);
});
if (reverse && penta_solve_method==penta_bicgstab) {
 // Keep the native optional initial-guess permutation. The current lower
 // solver rejects BiCGStab with the same supported scope as JAX.
 Real3D old_input(Kokkos::view_alloc(Kokkos::WithoutInitializing,"xm_wpxp_solve::old_input"),ngrdcol,2*nzm-1,nrhs);
 Kokkos::deep_copy(old_input,old_solution);
 Kokkos::parallel_for("xm_wpxp_solve::old_order",ngrdcol*(2*nzm-1),KOKKOS_LAMBDA(int cell) {
  int i=cell%ngrdcol,k=cell/ngrdcol;
  for(int n=0;n<nrhs;++n)old_solution(i,k,n)=old_input(i,2*nzm-2-k,n);
 });
}
Real3D rhs_solver=to_execution_space(rhs),solution_solver(Kokkos::view_alloc(Kokkos::WithoutInitializing,"xm_wpxp_solve::solution"),ngrdcol,2*nzm-1,nrhs);
Real3D old_solver=to_execution_space(old_solution);
Real1D condition;
if(rcond)condition=Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing,"xm_wpxp_solve::rcond"),ngrdcol);
band_solve("xm_wpxp",penta_solve_method,ngrdcol,nsup,nsub,2*nzm-1,nrhs,
           l_implemented,lhs_solver,rhs_solver,err_info,solution_solver,
           &old_solver,rcond?&condition:nullptr);
from_execution_space(rhs,rhs_solver);
if(rcond)from_execution_space(*rcond,condition);
Kokkos::parallel_for("xm_wpxp_solve::restore",ngrdcol*(2*nzm-1),KOKKOS_LAMBDA(int cell) {
 int i=cell%ngrdcol,k=cell/ngrdcol;
 for(int b=0;b<nsup+nsub+1;++b)lhs(b,i,k)=lhs_solver(i,b,k);
 for(int n=0;n<nrhs;++n)solution(i,k,n)=solution_solver(i,reverse?2*nzm-2-k:k,n);
});
if(f2k::clubb_at_least_debug_level_api(0)) {
 if(any_equal(err_info.err_code,clubb_fatal_error)) {
  std::cerr<<err_info.err_header_global<<"\\nError in xm_wpxp_solve\\n";
  return;
 }
}
'''

def correction_body(e):
 ns=bodies['xm_correction_wpxp_cl'].copy()
 first=next(i for i,n in enumerate(ns) if n.kind=='for')
 ns[first]=Node('raw','''// Native leaves untouched per-column mask entries undefined. Define them
// false as in JAX; no Native BFB claim applies to those undefined entries.
// Reduce integer column hits and keep each column's vertical scan ordered.
int clipping_columns=0;
Kokkos::parallel_reduce("xm_correction_wpxp_cl::clipping_columns",ngrdcol,
 KOKKOS_LAMBDA(int column,int& count) {
  bool needed=false;
  for(int level=0;level<nzm;++level) {
   if(Kokkos::abs(wpxp_chnge(column,level))>eps)needed=true;
  }
  l_clipping_needed(column)=needed;
  if(needed)++count;
 },clipping_columns);
l_any_clipping_needed=clipping_columns!=0;''')
 second=next(i for i,n in enumerate(ns) if n.kind=='for')
 ns.insert(second,Node('raw','''// The source also leaves nonclipping-column tacl statistics unwritten.
// Define their contribution as zero, matching JAX; physical updates below
// retain the Native masked formula and early-return placement.
Kokkos::deep_copy(xm_tndcy_wpxp_cl,zero);'''))
 return e.emit(ns)

support='''// Advance CLUBB means and turbulent fluxes with the native coupled equations.
//
// advance_xm_wpxp updates total water, liquid-water potential temperature,
// passive scalars, and optionally horizontal winds and their turbulent fluxes.
// Inputs/outputs use column-first execution-space Views. Private matrix arrays
// retain Native bands-first storage; the existing solver boundary transposes it.
// Stats retain the Fortran backend. Sponge damping is excluded with JAX scope
// and requires upstream rejection before this routine; see INTEGRATION.md.
// References: Golaz et al. (2002), JAS 59, 3540–3551, Eqs.16–17;
// https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wpxp_eqns
#include "advance_xm_wpxp_module.hpp"
#include "constants_clubb.hpp"
#include "parameter_indices.hpp"
#include "model_flags.hpp"
#include "clip_explicit.hpp" // Pending existing Native constant owner; see INTEGRATION.md.
#include "Skx_module.hpp"
#include "matrix_solver_wrapper.hpp"
#include "error_code.hpp"
#include "f2k/procedure_bridge.hpp"
#include "f2k_dependencies.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <algorithm>
namespace clubb_kokkos {
namespace {
'''
for k,v in module_constants.items():support+='constexpr int '+k+'='+v+';\n'
# These adapters implement unavoidable Native array reductions and host I/O;
# they contain no field-specific physics or diagnostic call lists.
support+='''template<class View> bool any_less(View values,real_t threshold) {
 // Express native ANY(values<threshold) across device cells for a host branch.
 int count=0;
 Kokkos::parallel_reduce("advance_xm_wpxp::any_less",values.extent(0)*values.extent(1),
 KOKKOS_LAMBDA(int cell,int& hits) {if(values(cell%values.extent(0),cell/values.extent(0))<threshold)++hits;},count);
 return count!=0;
}
template<class View>std::string dump_values(View values) {
 // Format a Native diagnostic array in first-index-fastest element order.
 auto host=to_host_space(values);std::ostringstream out;
 out<<std::setprecision(std::numeric_limits<real_t>::max_digits10);
 for(std::size_t n=0;n<host.size();++n)out<<host.data()[n]<<' ';
 return out.str();
}
template<class View>std::string dump_header(View values,int column) {
 // Convert one fixed-width Native character diagnostic to a host string.
 auto host=to_host_space(values);std::string result;
 for(std::size_t j=0;j<host.extent(0);++j)result+=host(j,column);
 return result;
}
'''
for n in names[1:]:support+=signature(n,True)+';\n'
support+='} // namespace\n\n'
for n in names:
 if n==names[1]:support+='namespace {\n'
 e=Emitter(n)
 support+=signature(n)+' {\n'
 support+='\n'.join('//'+t.text for t in headers[n])+'\n'
 for key,a in metadata[n].items():
  if a.name in args[n]:continue
  if a.dims:
   typ=scalar_type(a)+'*'*len(a.dims);typ='std::uint8_t'+'*'*len(a.dims) if a.kind=='logical' else typ
   support+='Kokkos::View<'+typ+',Kokkos::LayoutLeft> '+a.name+'(Kokkos::view_alloc(Kokkos::WithoutInitializing,"'+n+'::'+a.name+'"),'+','.join(e.expr(dim) for dim in a.dims)+');'+(' // '+variable_docs[n][key] if key in variable_docs[n] else '')+'\n'
  else:support+=('constexpr ' if key in initializers[n] else '')+scalar_type(a)+' '+a.name+('='+e.expr(initializers[n][key]) if key in initializers[n] else '')+';'+(' // '+variable_docs[n][key] if key in variable_docs[n] else '')+'\n'
 if n=='xm_wpxp_solve':body=solver_body(e)
 elif n=='xm_correction_wpxp_cl':body=correction_body(e)
 else:body=e.emit(bodies[n])+'\n'
 for key,alias in e.aliases.items():
  if key.startswith('err_info.err_header'):continue
  if re.search(r'\b'+alias+r'\b',re.sub(r'//[^\n]*','',body)):
   support+='const auto '+alias+'='+key+';\n'
 support+=body
 support+='}\n\n'
support+='} // namespace\n} // namespace clubb_kokkos\n'
head='''// Public Native-order mean/flux timestep interface; see the source description.
#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"
#include "parameters_tunable.hpp"
#include "pdf_parameter_module.hpp"
#include "err_info_type_module.hpp"
#include "f2k/stats_bridge.hpp"
namespace clubb_kokkos {
'''+signature(names[0])+';\n} // namespace clubb_kokkos\n'
(folder/'advance_xm_wpxp_module.cpp').write_text(support)
(folder/'advance_xm_wpxp_module.hpp').write_text(head)
