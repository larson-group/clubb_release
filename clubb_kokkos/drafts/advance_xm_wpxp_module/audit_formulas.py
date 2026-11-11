"""Compare explicit Native assignment operation trees without executing physics."""
from pathlib import Path
import re,sys,ast,json,hashlib
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import statements,declaration
source=(root/'src/CLUBB_core/advance_xm_wpxp_module.F90').read_text()
target=(folder/'advance_xm_wpxp_module.cpp').read_text()
# Reuse contract-audit parsing/canonical integer index representation only.
# Its top level also verifies the current final source/call/format artifacts.
import audit_source as contract
cpp=contract.cpp
counts={};skipped={}

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

for name in contract.native_names:
 if name=='xm_wpxp_solve':skipped[name]='Storage/solver API boundary audited separately';continue
 native=contract.nb[name];actual=contract.cb[name]
 meta={}
 for s in statements(native):
  if '::' in s:meta.update(declaration(s))
 arrays={a.name for a in meta.values() if a.dims}
 arrays.update(re.findall(r'\b(gr%\w+|pdf_implicit_coefs_terms%\w+|err_info%err_code)\s*\(',native))
 arrays={a.replace('%','.') for a in arrays}
 aliases={m[1]:m[2] for m in re.finditer(r'const auto\s+(\w+)\s*=\s*([\w.]+)\s*;',actual)}
 pattern=r'(?<![\w.])('+'|'.join(re.escape(a) for a in sorted(arrays,key=len,reverse=True))+r')\s*\('
 def expr(s,is_cpp=False,as_key=False):
  s=re.sub(r'\s+',' ',s.strip()).replace('%','.')
  for alias,original in aliases.items():s=re.sub(r'\b'+alias+r'\b',original,s)
  if not is_cpp:s=square_expand(s)
  s=re.sub(r'_core_rknd\b','',s).replace('Kokkos::','')
  s=re.sub(r'\breal_t\(([^()]*)\)',r'(\1)',s)
  for original,new in [('.true.','True'),('.false.','False'),('true','True'),('false','False')]:s=s.replace(original,new)
  for original,new in [('.and.',' and '),('.or.',' or '),('.not.',' not '),('&&',' and '),('||',' or ')]:s=s.replace(original,new)
  refs=[]
  while m:=re.search(pattern,s):
   raw,end=parens(s,m.end());indices=contract.split_top(raw)
   if any(':' in i for i in indices):raise ValueError('whole/section assignment')
   normalized=[]
   for index in indices:
    if index.startswith(('iC','ic_','ialtitude','iwpxp')):normalized.append(index)
    else:normalized.append(contract.polynomial('('+index+')+1' if is_cpp else index))
   ref=(m[1].lower(),tuple(normalized));token='ARRAY_'+str(len(refs));refs.append(ref)
   s=s[:m.start()]+token+s[end:]
  if as_key:
   assert len(refs)==1,(name,s,refs)
   return refs[0]
  node=ast.parse(s,mode='eval').body
  def norm(n):
   if isinstance(n,ast.Name):return refs[int(n.id[6:])] if n.id.startswith('ARRAY_') else ('name',n.id.lower())
   if isinstance(n,ast.Constant):return ('constant',str(n.value) if isinstance(n.value,str) else float(n.value).hex())
   if isinstance(n,ast.UnaryOp):return (type(n.op).__name__,norm(n.operand))
   if isinstance(n,ast.BinOp):return (type(n.op).__name__,norm(n.left),norm(n.right))
   if isinstance(n,ast.BoolOp):return (type(n.op).__name__,tuple(norm(v) for v in n.values))
   if isinstance(n,ast.Call):return ('call',n.func.id.lower(),tuple(norm(a) for a in n.args))
   raise ValueError(ast.dump(n))
  return norm(node)

 expected={};ignored=0
 for st in statements(native):
  m=re.match(r'^([\w%]+)(\(.*?\))?\s*=(?!=)\s*(.*)$',st)
  if not m:continue
  var,index,rhs=m.groups();var=var.replace('%','.')
  if var not in arrays and var.lower() not in meta:continue
  if var in {'i','k','j','n','sclr','l_any_clipping_needed','l_clipping_needed'}:continue
  if var in arrays and (not index or ':' in index):continue
  if 'sponge_damp' in rhs:continue
  if re.search(r'\b(?:max|min|present)\s*\(',rhs,re.I):ignored+=1;continue
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
  if '?' in rhs:continue
  # Ignore declarations and loop-coordinate/private-parameter initializations;
  # explicit body variables have one additional token only when declared.
  prefix=actual[max(0,m.start()-30):m.start()]
  if re.search(r'\b(?:constexpr|const)\s+(?:real_t|bool|int)\s*$',prefix):continue
  try:observed[key].append(expr(rhs,True))
  except (ValueError,SyntaxError):continue
 for key,values in expected.items():
  assert observed[key]==values,(name,key,values,observed[key])
 counts[name]=sum(len(v) for v in expected.values());skipped[name]=ignored
record={'explicit_assignment_operation_trees':counts,'total':sum(counts.values()),
        'manual_minmax_or_storage_cases':skipped,
        'native_sha256':hashlib.sha256(source.encode()).hexdigest(),
        'target_sha256':hashlib.sha256(target.encode()).hexdigest(),
        'scope':'Source operation structure only; no compilation or numeric execution'}
(folder/'formula_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
