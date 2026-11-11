"""Compare explicit Native assignment operation trees without executing physics."""
from pathlib import Path
import re,sys,ast,json,hashlib
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import statements,declaration
source=(root/'src/CLUBB_core/advance_xp2_xpyp_module.F90').read_text()
target=(folder/'advance_xp2_xpyp_module.cpp').read_text()
# Reuse contract-audit parsing/canonical integer index representation only.
# Its top level also verifies the current final source/call/format artifacts.
import audit_source as contract
cpp=contract.cpp
counts={};skipped={};branch_counts={};whole_counts={}

def parens(s,start):return contract.balanced_at(s,start)

def square_expand(s):
 while m:=re.search(r'\*\*\s*2\b',s):
  end=m.start();start=end-1
  while s[start].isspace():start-=1
  if s[start]==')':
   depth=1;start-=1
   while depth:depth+=(s[start]==')')-(s[start]=='(');start-=1
  while start>=0 and (s[start].isalnum() or s[start] in '_.%'):start-=1
  start+=1;base=s[start:end];s=s[:start]+'('+base+'*'+base+')'+s[m.end():]
 return s

def ternaries(text):
    """Translate fully parenthesized C++ conditional expressions to Python AST syntax."""
    result='';i=0
    while i<len(text):
        if text[i]=='(':
            depth=1;end=i+1
            while depth:
                if text[end]=='(':depth+=1
                elif text[end]==')':depth-=1
                end+=1
            result+='('+ternaries(text[i+1:end-1])+')';i=end
        else:result+=text[i];i+=1
    depth=0;question=None;colon=None
    for i,ch in enumerate(result):
        if ch=='(':depth+=1
        elif ch==')':depth-=1
        elif depth==0 and ch=='?':question=i
        elif depth==0 and ch==':' and question is not None:colon=i;break
    if question is not None:
        if colon is None:raise ValueError(('ternary colon',result))
        result='('+result[question+1:colon]+' if '+result[:question]+' else '+result[colon+1:]+')'
    return result



for name in contract.native_names:
 if name=='xp2_xpyp_solve':skipped[name]='Storage/solver API boundary audited separately';continue
 native=contract.nb[name];actual=contract.cb[name]
 meta={}
 for s in statements(native):
  if '::' in s:meta.update(declaration(s))
 arrays={a.name for a in meta.values() if a.dims}
 arrays.update(re.findall(r'\b(gr%\w+|pdf_implicit_coefs_terms%\w+|pdf_params%\w+|err_info%err_code)\s*\(',native))
 arrays={a.replace('%','.') for a in arrays}
 aliases={m[1]:m[2] for m in re.finditer(r'const auto\s+(\w+)\s*=\s*([\w.]+)\s*;',actual)}
 pattern=r'(?<![\w.])('+'|'.join(re.escape(a) for a in sorted(arrays,key=len,reverse=True))+r')\s*\('
 def expr(s,is_cpp=False,as_key=False):
  s=re.sub(r'\s+',' ',s.strip()).replace('%','.').replace('/=','!=')
  s=s.replace('stats.flags().l_sample','stats.l_sample')
  s=re.sub(r'any\(\s*err_info.err_code\s*==\s*clubb_fatal_error\s*\)','any_equal(err_info.err_code,clubb_fatal_error)',s,flags=re.I)
  for alias,original in aliases.items():s=re.sub(r'\b'+alias+r'\b',original,s)
  if is_cpp:s=re.sub(r'\bcolumn\b','(i-1)',s)
  if not is_cpp:s=square_expand(s)
  s=re.sub(r'_core_rknd\b','',s).replace('Kokkos::','').replace('f2k::','').replace('copysign(', 'sign(')
  s=re.sub(r'\breal_t\(([^()]*)\)',r'(\1)',s)
  for original,new in [('.true.','True'),('.false.','False'),('true','True'),('false','False')]:s=s.replace(original,new)
  for original,new in [('.and.',' and '),('.or.',' or '),('.not.',' not '),('&&',' and '),('||',' or ')]:s=s.replace(original,new)
  refs=[]
  while m:=re.search(pattern,s):
   raw,end=parens(s,m.end());indices=contract.split_top(raw)
   if any(':' in i for i in indices):raise ValueError('whole/section assignment')
   normalized=[]
   for index in indices:
    if index.strip() in contract.parameter_indices:normalized.append(index)
    else:normalized.append(contract.polynomial('('+index+')+1' if is_cpp else index))
   ref=(m[1].lower(),tuple(normalized));token='ARRAY_'+str(len(refs));refs.append(ref)
   s=s[:m.start()]+token+s[end:]
  if as_key:
   assert len(refs)==1,(name,s,refs)
   return refs[0]
  s=re.sub(r'!(?!=)',' not ',s)
  s=ternaries(s.strip())
  node=ast.parse(s,mode='eval').body
  def norm(n):
   if isinstance(n,ast.Name):return refs[int(n.id[6:])] if n.id.startswith('ARRAY_') else ('name',n.id.lower())
   if isinstance(n,ast.Constant):return ('constant',str(n.value) if isinstance(n.value,str) else float(n.value).hex())
   if isinstance(n,ast.UnaryOp):return (type(n.op).__name__,norm(n.operand))
   if isinstance(n,ast.BinOp):return (type(n.op).__name__,norm(n.left),norm(n.right))
   if isinstance(n,ast.BoolOp):return (type(n.op).__name__,tuple(norm(v) for v in n.values))
   if isinstance(n,ast.Attribute):return ('name',ast.unparse(n).lower())
   if isinstance(n,ast.Compare):return ('compare',norm(n.left),tuple(type(op).__name__ for op in n.ops),tuple(norm(v) for v in n.comparators))
   if isinstance(n,ast.IfExp):
    test=n.test
    if isinstance(test,ast.Compare) and len(test.ops)==1 and norm(n.body)==norm(test.left) and norm(n.orelse)==norm(test.comparators[0]):return ('call','max' if isinstance(test.ops[0],ast.Gt) else 'min',(norm(n.body),norm(n.orelse)))
    raise ValueError(ast.dump(n))
   if isinstance(n,ast.Call):
    if isinstance(n.func,ast.Attribute) and n.func.attr=='size':return ('call','len_trim',(norm(n.func.value),))
    return ('call',n.func.id.lower(),tuple(norm(a) for a in n.args))
   raise ValueError(ast.dump(n))
  return norm(node)

 expected={};ignored=0
 for st in statements(native):
  m=re.match(r'^([\w%]+)(\(.*?\))?\s*=(?!=)\s*(.*)$',st)
  if not m:continue
  var,index,rhs=m.groups();var=var.replace('%','.')
  if var not in arrays and var.lower() not in meta:continue
  if var in {'i','k','j','n','sclr','l_any_clipping_needed','l_clipping_needed','l_single_solve_possible'}:continue
  if var in arrays and (not index or ':' in index):continue
  if 'sponge_damp' in rhs:continue
  if re.search(r'\bpresent\s*\(',rhs,re.I):ignored+=1;continue
  key=expr(var+index,as_key=True) if index else ('scalar',var.lower())
  expected.setdefault(key,[]).append(expr(rhs))
 observed={key:[] for key in expected}
 variants={a.name for a in meta.values()}|set(aliases)
 for m in re.finditer(r'(?<![\w.])('+'|'.join(re.escape(a) for a in sorted(variants,key=len,reverse=True))+r')\b',actual):
  var=m[1];pos=m.end();index=''
  while pos<len(actual) and actual[pos].isspace():pos+=1
  if actual[pos:pos+1]=='(':
   raw,end=parens(actual,pos+1);index='('+raw+')';pos=end
  while pos<len(actual) and actual[pos].isspace():pos+=1
  if actual[pos:pos+1]!='=' or actual[pos:pos+2]=='==':continue
  end=actual.index(';',pos)
  rhs=actual[pos+1:end].strip()
  try:key=expr(var+index,True,True) if index else ('scalar',var.lower())
  except (ValueError,SyntaxError):continue
  if key not in expected:continue
  # Ternaries are normalized below without discarding MIN/MAX operands.
  # Ignore declarations and loop-coordinate/private-parameter initializations;
  # explicit body variables have one additional token only when declared.
  prefix=actual[max(0,m.start()-30):m.start()]
  if var not in {'C2rt','C2thl','C2rtthl'} and re.search(r'\b(?:constexpr|const)\s+(?:real_t|bool|int)\s*$',prefix):continue
  try:observed[key].append(expr(rhs,True))
  except (ValueError,SyntaxError):continue
 for key,values in expected.items():
  assert observed[key]==values,(name,key,values,observed[key])
 # Independently expand whole-array expression assignments into logical slice
 # coordinates, then compare their operation trees and every affine index.
 ranks={a.name:len(a.dims) for a in meta.values() if a.dims}
 for component in arrays:
  if component in ranks:continue
  if m:=re.search(re.escape(component.replace('.', '%'))+r'\s*\(',native):
   raw,end=parens(native,m.end());ranks[component]=len(contract.split_top(raw))
 def expand_whole(text):
  text=text.replace('%','.');protected=[]
  pattern_all=r'(?<![\w.])('+'|'.join(re.escape(a) for a in sorted(ranks,key=len,reverse=True))+r')\b'
  out='';pos=0
  while m:=re.search(pattern_all,text[pos:]):
   start=pos+m.start();end=pos+m.end();var=m[1];out+=text[pos:start];cursor=end
   while cursor<len(text) and text[cursor].isspace():cursor+=1
   if text[cursor:cursor+1]=='(':
    raw,end=parens(text,cursor+1);indices=contract.split_top(raw);coordinate=0;expanded=[]
    for index in indices:
     if ':' in index:
      lo,*parts=index.split(':');lo=lo.strip() or '1';step=parts[1].strip() if len(parts)>1 else '1'
      expanded.append('('+lo+'+slice_'+str(coordinate)+'*('+step+'))');coordinate+=1
     else:expanded.append(index)
   else:expanded=['1+slice_'+str(j) for j in range(ranks[var])]
   out+=var+'('+','.join(expanded)+')';pos=end
  return out+text[pos:]
 whole_expected={};whole_observed={}
 for st in statements(native):
  m=re.match(r'^([\w%]+)(\(.*?\))?\s*=(?!=)\s*(.*)$',st)
  if not m:continue
  var,index,rhs=m.groups();var=var.replace('%','.')
  if var not in ranks or index is not None and ':' not in index:continue
  if re.match(r'\s*(?:zt2zm_api|zm2zt_api|sponge_damp_xp2)\s*\(',rhs):continue
  key=expr(expand_whole(var+(index or '')),as_key=True)
  whole_expected.setdefault(key,[]).append(expr(expand_whole(rhs)))
 whole_observed={key:[] for key in whole_expected}
 variants=set(ranks)|set(aliases)
 for m in re.finditer(r'(?<![\w.])('+'|'.join(re.escape(a) for a in sorted(variants,key=len,reverse=True))+r')\s*\(',actual):
  raw,pos=parens(actual,m.end());tail=pos
  while tail<len(actual) and actual[tail].isspace():tail+=1
  if actual[tail:tail+1]!='=' or actual[tail:tail+2]=='==':continue
  end=actual.index(';',tail);value=actual[tail+1:end]
  try:key=expr(m[1]+'('+raw+')',True,True)
  except (ValueError,SyntaxError):continue
  if key not in whole_expected:continue
  whole_observed[key].append(expr(value,True))
 for key,values in whole_expected.items():
  assert whole_observed[key]==values,(name,'whole-array assignment',key,values,whole_observed[key])
 whole_counts[name]=sum(len(v) for v in whole_expected.values())
 # Compare every actual Native IF predicate, including one-line diagnostic returns.
 # The solver's storage-boundary branches are reviewed separately above.
 native_conditions=[];sponge=0
 for st in statements(native):
  if re.match(r'if\s*\(',st,re.I) and re.search(r'\bthen\s*$',st,re.I):
   if 'sponge_damp_settings%' in st:sponge=1;continue
   if sponge:sponge+=1
  if re.match(r'end\s*if',st,re.I) and sponge:sponge-=1;continue
  if sponge:continue
  if m:=re.match(r'(?:if|elseif|else\s+if)\s*\(',st,re.I):
   raw,end=parens(st,m.end());native_conditions.append(expr(raw))
 actual_conditions=[]
 for m in re.finditer(r'\bif\s*\(',actual):
  raw,end=parens(actual,m.end());actual_conditions.append(expr(raw,True))
 assert actual_conditions==native_conditions,(name,'condition order',native_conditions,actual_conditions)
 # Native SELECT CASE identifiers remain literal model-category names.
 native_cases=[re.sub(r'\s+','',m[1]).lower() for m in re.finditer(r'^\s*case\s*\(([^)]*)\)',native,re.M|re.I)]
 actual_cases=[re.sub(r'\s+','',m[1]).lower() for m in re.finditer(r'\bcase\s+([^:]+):',actual)]
 assert native_cases==actual_cases,(name,'case order',native_cases,actual_cases)
 counts[name]=sum(len(v) for v in expected.values());skipped[name]=ignored
 branch_counts[name]=len(native_conditions)

record={'explicit_assignment_operation_trees':counts,'total':sum(counts.values()),'whole_array_operation_trees':whole_counts,'whole_array_total':sum(whole_counts.values()),'ordered_if_predicates':branch_counts,'total_if_predicates':sum(branch_counts.values()),
        'manual_minmax_or_storage_cases':skipped,
        'native_sha256':hashlib.sha256(source.encode()).hexdigest(),
        'target_sha256':hashlib.sha256(target.encode()).hexdigest(),
        'scope':'Source operation structure only; no compilation or numeric execution'}
(folder/'formula_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
