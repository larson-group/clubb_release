"""Audit third-moment source contracts and expressions without running physics."""
from pathlib import Path
import ast,re,sys,json,hashlib
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import read_procedure,statements,declaration
from generate_derived_types import split_top
path=root/'src/CLUBB_core/advance_xp3_module.F90';source=path.read_text()
target=(folder/'advance_xp3_module.cpp').read_text();header=(folder/'advance_xp3_module.hpp').read_text()
names=['advance_xp3','diagnose_xp3','advance_xp3_simplified','term_tp_rhs','term_ac_rhs']
clean=lambda s:re.sub(r'//[^\n]*','',s)
cpp=clean(target);nb={};cb={};metadata={};locations=[]

def bracket(s,start):
 depth,pos=1,start
 while depth:depth+=(s[pos]=='(')-(s[pos]==')');pos+=1
 return s[start:pos-1],pos

def brace(s,start):
 depth,pos=1,start
 while depth:depth+=(s[pos]=='{')-(s[pos]=='}');pos+=1
 return s[start:pos-1]

for name in names:
 nb[name]=re.search(r'^  (?:subroutine|function) '+name+r'\(.*?^  end (?:subroutine|function) '+name+r'\b',source,re.M|re.S).group()
 match=list(re.finditer(r'(?:void|real_t)\s+'+name+r'\(([^()]*)\)\s*\{',cpp))
 assert len(match)==1,(name,len(match))
 m=match[0];locations.append(m.start());cb[name]=brace(cpp,m.end())
 procedure=read_procedure(path,name)
 args=[part.split()[-1].lstrip('&*') for part in split_top(m[1])]
 assert args==[a.name for a in procedure.arguments]
 assert not any(a.optional for a in procedure.arguments)
 if name in names[:2]:
  h=re.search(r'void\s+'+name+r'\(([^()]*)\)',clean(header))
  assert [part.split()[-1].lstrip('&*') for part in split_top(h[1])]==args
 metadata[name]={}
 for st in statements(nb[name]):
  if '::' in st:metadata[name].update(declaration(st))
assert locations==sorted(locations)
all_calls=set(names)|{'zm2zt_api','zt2zm_api','ddzm','Skx_func','xp3_LG_2005_ansatz','stats_update','stats_begin_budget','stats_finalize_budget'}

def calls(text):
 result=[]
 for m in re.finditer(r'(?<!\w)(?:f2k::)?('+'|'.join(sorted(all_calls,key=len,reverse=True))+r')\s*\(',text):
  raw,end=bracket(text,m.end());result.append((m[1],split_top(raw)))
 return result


def square_expand(s):
 while m:=re.search(r'\*\*\s*2\b',s):
  end=m.start();start=end-1
  while s[start].isspace():start-=1
  if s[start]==')':
   d=1;start-=1
   while d:d+=(s[start]==')')-(s[start]=='(');start-=1
  while start>=0 and (s[start].isalnum() or s[start] in '_.%'):start-=1
  start+=1;b=s[start:end];s=s[:start]+'('+b+'*'+b+')'+s[m.end():]
 return s


def polynomial(s):
 node=ast.parse(s,mode='eval').body
 def add(a,b,scale=1):
  out=dict(a)
  for k,v in b.items():out[k]=out.get(k,0)+scale*v
  return {k:v for k,v in out.items() if v}
 def norm(n):
  if isinstance(n,ast.Name):return {n.id:1}
  if isinstance(n,ast.Constant):return {'1':n.value}
  if isinstance(n,ast.UnaryOp):return {k:(-v if isinstance(n.op,ast.USub) else v) for k,v in norm(n.operand).items()}
  if isinstance(n,ast.BinOp):
   a,b=norm(n.left),norm(n.right)
   if isinstance(n.op,ast.Add):return add(a,b)
   if isinstance(n.op,ast.Sub):return add(a,b,-1)
   if isinstance(n.op,ast.Mult):
    if set(a)<={'1'}:return {k:a.get('1',0)*v for k,v in b.items()}
    if set(b)<={'1'}:return {k:b.get('1',0)*v for k,v in a.items()}
  raise ValueError(s)
 return tuple(sorted(norm(node).items()))


def expression(raw,name,is_cpp=False,whole=False):
 s=re.sub(r'\s+',' ',raw.strip()).replace('%','.')
 s=s.replace('sclr_tol_host','sclr_tol').replace('invrs_dzt(', 'gr.invrs_dzt(') if is_cpp else s
 if not is_cpp:s=square_expand(s)
 s=re.sub(r'_core_rknd\b','',s).replace('Kokkos::','')
 s=re.sub(r'\breal_t\(([^()]*)\)',r'(\1)',s)
 if s=='sclrp3_field':s='subview(sclrp3,ALL,ALL,sclr)'
 arrays={a.name for a in metadata[name].values() if a.dims}
 arrays.add('gr.invrs_dzt')
 refs=[]
 def make_ref(base,indices):
  if all(index==':' for index in indices):return base
  key=[]
  for index in indices:
   if index==':':key.append(':');continue
   if index.startswith(('ixp3_','iC')):key.append(index);continue
   if not is_cpp:index=re.sub(r'\b(i|k|kp1|sclr)\b',r'(\1+1)',index);index='('+index+')-1'
   key.append(polynomial(index))
  token='ARRAY_'+str(len(refs));refs.append((base.lower(),tuple(key)));return token
 while m:=re.search(r'\bsubview\(',s):
  body,end=bracket(s,m.end());parts=split_top(body)
  replacement=make_ref(parts[0],[(':' if p=='ALL' else p) for p in parts[1:]])
  s=s[:m.start()]+replacement+s[end:]
 pattern=r'(?<![\w.])('+'|'.join(re.escape(a) for a in sorted(arrays,key=len,reverse=True))+r')\s*\('
 while m:=re.search(pattern,s):
  body,end=bracket(s,m.end());indices=split_top(body)
  if whole and is_cpp and all(i in {'i','k'} for i in indices):replacement=m[1]
  else:replacement=make_ref(m[1],indices)
  s=s[:m.start()]+replacement+s[end:]
 node=ast.parse(s,mode='eval').body
 def norm(n):
  if isinstance(n,ast.Name):return refs[int(n.id[6:])] if n.id.startswith('ARRAY_') else ('name',n.id.lower())
  if isinstance(n,ast.Constant):return ('constant',n.value if isinstance(n.value,str) else float(n.value).hex())
  if isinstance(n,ast.UnaryOp):return (type(n.op).__name__,norm(n.operand))
  if isinstance(n,ast.BinOp):return (type(n.op).__name__,norm(n.left),norm(n.right))
  if isinstance(n,ast.Call):return ('call',n.func.id.lower(),tuple(norm(a) for a in n.args))
  raise ValueError(ast.dump(n))
 return norm(node)

counts={};checked_args=0
for name in names:
 ss=statements(nb[name]);lastdecl=max(i for i,s in enumerate(ss) if '::' in s)
 original=calls('\n'.join(ss[lastdecl+1:-1]));actual=calls(cb[name])
 assert [n for n,_ in original]==[n for n,_ in actual],(name,[n for n,_ in original],[n for n,_ in actual])
 counts[name]=len(actual)
 bindings={}
 for m in re.finditer(r'\b(stats_bt_before|stats_dp|stats_bt_after)\(i,\s*k\)\s*=\s*([^;]+);',cb[name]):bindings[m[1]]=m[2]
 for (native_name,expected),(cpp_name,got) in zip(original,actual):
  assert len(expected)==len(got),(name,native_name)
  for j,(e,a) in enumerate(zip(expected,got)):
   whole=cpp_name.startswith('stats_') and j==1 and a in bindings
   if whole:a=bindings[a]
   try:assert expression(e,name)==expression(a,name,True,whole)
   except Exception:raise AssertionError((name,native_name,j,e,a))
   checked_args+=1

# Exact scalar helper association, private update and coefficient formulas.
formula_checks=0
for name,result in [('term_tp_rhs','term_tp'),('term_ac_rhs','term_ac')]:
 statements_=statements(nb[name]);rhs=next(s.split('=',1)[1] for s in statements_ if re.match(result+r'\s*=',s))
 actual=re.search(r'\breturn\s+([^;]+);',cb[name])[1]
 assert expression(rhs,name)==expression(actual,name,True)
 formula_checks+=1
for name in names[:2]:
 rhs=next(s.split('=',1)[1] for s in statements(nb[name]) if re.match(r'xp3_coef_fnc\(i,k\)\s*=',s))
 actual=re.search(r'xp3_coef_fnc\(i,\s*k\)\s*=\s*([^;]+);',cb[name])[1]
 assert expression(rhs,name)==expression(actual,name,True)
 formula_checks+=1
 rhs=next(s.split('=',1)[1] for s in statements(nb[name]) if s.startswith('brunt_vaisala_freq_sqd_zt(:,:)'))
 stability=re.search(r'const real_t stability\s*=\s*([^;]+);',cb[name])[1]
 assert re.search(r'brunt_vaisala_freq_sqd_zt\(i,\s*k\)\s*=\s*stability > zero \? stability : zero',cb[name])
 assert expression(rhs,name)==expression('max('+stability+',zero)',name,True,True)
 formula_checks+=1
rhs=next(s.split('=',1)[1] for s in statements(nb['advance_xp3_simplified']) if s.startswith('xp3(i,k) = ( ( xp3'))
actual=re.search(r'xp3\(i,\s*k\)\s*=\s*(\(\(xp3[^;]+);',cb['advance_xp3_simplified'])[1]
assert expression(rhs,'advance_xp3_simplified')==expression(actual,'advance_xp3_simplified',True)
formula_checks+=1
rhs=next(s.split('=',1)[1] for s in statements(nb['advance_xp3_simplified']) if s.startswith('xp3(i,k) = min('))
actual=re.search(r'xp3\(i,\s*k\)\s*=\s*(\(reciprocal_tau[^;]+);',cb['advance_xp3_simplified'])[1]
actual=re.sub(r'\(reciprocal_tau < tau_max_zt\(i,\s*k\) \? reciprocal_tau : tau_max_zt\(i,\s*k\)\)',
              'min((one/invrs_tau_zt(i,k)),tau_max_zt(i,k))',actual)
assert expression(rhs,'advance_xp3_simplified')==expression(actual,'advance_xp3_simplified',True)
formula_checks+=1
rhs=next(s.split('=',1)[1] for s in statements(nb['advance_xp3_simplified']) if s.startswith('xp3(i,k) = one_half'))
actual=re.search(r'xp3\(i,\s*k\)\s*=\s*(one_half[^;]+);',cb['advance_xp3_simplified'])[1]
assert expression(rhs,'advance_xp3_simplified')==expression(actual,'advance_xp3_simplified',True)
formula_checks+=1

# Branchs, special minima, source boundary and ownership are separately visible.
assert 'constexpr bool l_predict_xp3 = false' in target
assert 'constexpr real_t C_xp3_dissipation = real_t(1.0)' in target
assert 'xp3(i, nzt - 1) = zero' in target
assert 'const int kp1 = k + 1 < nzt - 1 ? k + 1 : nzt - 1' in target
assert 'reciprocal_tau < tau_max_zt(i, k) ? reciprocal_tau : tau_max_zt(i, k)' in target
assert 'if (l_predict_xp3 && l_lmm_stepping)' in target
assert 'if (iiPDF_type == iiPDF_ADG1)' in target
assert target.count('from_execution_space(sclrp3_section, sclrp3_field)')==2
assert target.count('to_execution_space(sclrp3_section, false)')==2
assert 'f2k::advance_xp3' not in target and 'f2k::diagnose_xp3' not in target
assert 'fma(' not in cpp
kernels=0
for m in re.finditer(r'KOKKOS_LAMBDA\([^{}]*?\)\s*\{',cpp):
 b=brace(cpp,m.end());kernels+=1
 assert not re.search(r'\bgr\.',b)
 assert not any(token in b for token in ['f2k::','to_host_space(','Kokkos::parallel_for('])
for filename in ['advance_xp3_module.cpp','advance_xp3_module.hpp']:
 text=(folder/filename).read_text();code=clean(text);code=re.sub(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'','',code)
 for left,right in [('(',')'),('[',']'),('{','}')]:
  depth=0
  for char in code:depth+=(char==left)-(char==right);assert depth>=0,(filename,left)
  assert depth==0,(filename,left)
 assert all(len(line)<=100 and line.rstrip()==line and '\t' not in line for line in text.splitlines())
record={
 'status':'Complete supported source draft; not compiled, integrated or numerically validated',
 'source':'src/CLUBB_core/advance_xp3_module.F90; all 1006 lines read',
 'routine_order':names,'argument_counts':{n:len(read_procedure(path,n).arguments) for n in names},
 'all_arguments_required_no_optionals':True,'call_counts_including_grid_and_scalar_functions':counts,
 'all_call_arguments_checked':checked_args,'explicit_formula_trees_checked':formula_checks,
 'kernel_bodies_checked':kernels,'no_owning_udt_host_reads_bridges_nested_launches_in_kernels':True,
 'manual_audit':['Full Native predictive-equation and discretization comments retained','ADG1/non-ADG1 scalar-flux minimum difference retained','Steady-state cap and buoyancy-frequency MAX use Native operand order','Native last-storage-level boundary retained independently of grid direction','Predictive/LMM snapshots, all source stats names/gates and mean/variance minima','All scratch is per-call; only existing array owner handles scalar/subview boundaries','No vertical recurrence/matrix solve in scope; no missing lower bridge/constant owner'],
 'missing_lower_dependencies':[],'runtime':'No builds, installs, tests or numerical claims; root owns CPU/CPU zero-tolerance model gate and occasional GPU checkpoints',
 'native_sha256':hashlib.sha256(source.encode()).hexdigest(),
 'sha256':{filename:hashlib.sha256((folder/filename).read_bytes()).hexdigest() for filename in ['advance_xp3_module.cpp','advance_xp3_module.hpp']}
}
(folder/'source_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
