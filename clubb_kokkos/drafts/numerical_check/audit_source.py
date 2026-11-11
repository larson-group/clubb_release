"""Independently audit the source draft without compiling or executing science."""
from pathlib import Path
import ast,hashlib,json,re,sys
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import statements
source=(root/'src/CLUBB_core/numerical_check.F90').read_text()
target=(folder/'numerical_check.cpp').read_text();header=(folder/'numerical_check.hpp').read_text()
inventory=json.loads((folder/'SOURCE_OUTLINE.json').read_text())['routines'];names=[r['name'] for r in inventory]
live={'parameterization_check','check_negative','check_nan_1d','check_nan_2d'}
new=[n for n in names if n not in live]
nb={n:re.search(r'^  (?:logical\s+)?(?:subroutine|function) '+n+r'\s*\(.*?^\s*end (?:subroutine|function) '+n+r'\b',source,re.M|re.S).group() for n in names}
# Preserve literals while stripping comments, including '%' and commas in labels.
literal=r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\''
def clean(text):return re.sub(literal+r'|//[^\n]*',lambda m:m[0] if not m[0].startswith('//') else ' '*len(m[0]),text)
def mask(text):return re.sub(literal,lambda m:' '*len(m[0]),clean(text))
def split_top(text,separator=','):
 result=[];quote='';depth=0;start=0;i=0
 while i<len(text):
  ch=text[i]
  if quote:
   if ch==quote and (i==0 or text[i-1]!='\\'):
    if i+1<len(text) and text[i+1]==quote:i+=1
    else:quote=''
  elif ch in "\"'":quote=ch
  elif ch in '([{':depth+=1
  elif ch in ')]}':depth-=1
  elif text[i:i+len(separator)]==separator and depth==0:result.append(text[start:i].strip());start=i+len(separator);i+=len(separator)-1
  i+=1
 result.append(text[start:].strip());return result

def balanced(text,start,left='(',right=')'):
 probe=mask(text);depth=1;pos=start
 while depth:depth+=(probe[pos]==left)-(probe[pos]==right);pos+=1
 return text[start:pos-1],pos

def definition(text,name):
 probe=mask(text);m=re.search(r'(?:(?:template <class View>\s*)?static\s+)?(?:void|bool|real_t)\s+'+name+r'\s*\([^{};]*\)\s*\{',probe)
 assert m,name
 body,end=balanced(text,m.end(),'{','}')
 signature=text[m.start():m.end()];return signature,body,text[m.start():end],m.start()

cb={};locations=[];counts={};preserved={}
for routine in inventory:
 n=routine['name'];container=header if n=='calculate_spurious_source' else target
 sig,body,full,offset=definition(container,n);cb[n]=clean(body)
 if n=='calculate_spurious_source':offset=target.index('// Native calculate_spurious_source')
 locations.append(offset)
 sig=clean(sig);raw=sig[sig.index('(')+1:sig.rindex(')')]
 actual=[a.split('=')[0].strip().split()[-1].lstrip('*&') for a in split_top(raw)]
 expected=[a['name'] for a in routine['arguments']]
 assert actual==expected,(n,actual,expected)
 assert all(not a['optional'] for a in routine['arguments']),n
 counts[n]=len(actual)
 if n in live:
  starting=(folder/'starting_live/numerical_check.cpp').read_text();old=definition(starting,n)[2]
  assert full==old,('live body changed',n)
  preserved[n]=hashlib.sha256(full.encode()).hexdigest()
 if routine['scope']=='public':
  h=clean(header);m=re.search(r'\b(?:void|bool|real_t)\s+'+n+r'\s*\(',h)
  args,_=balanced(h,m.end())
  assert [a.split('=')[0].strip().split()[-1].lstrip('*&') for a in split_top(args)]==expected,n
assert locations==sorted(locations),'Native definition order'
old_header=(folder/'starting_live/numerical_check.hpp').read_text()
old_public=old_header[old_header.index('void parameterization_check('):old_header.index('\n\n} // namespace')]
assert old_public in header

param_indices=set(re.findall(r'inline constexpr int (i\w+)',(root/'clubb_kokkos/src/CLUBB_core/parameter_indices.hpp').read_text()))

def polynomial(text):
 tree=ast.parse(text,mode='eval').body
 def norm(n):
  if isinstance(n,ast.Constant):return {'1':n.value}
  if isinstance(n,ast.Name):return {n.id:1}
  if isinstance(n,ast.UnaryOp):return {k:(-v if isinstance(n.op,ast.USub) else v) for k,v in norm(n.operand).items()}
  if isinstance(n,ast.BinOp):
   a,b=norm(n.left),norm(n.right)
   if isinstance(n.op,(ast.Add,ast.Sub)):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,0)+(-v if isinstance(n.op,ast.Sub) else v)
    return {k:v for k,v in out.items() if v}
   if isinstance(n.op,ast.Mult):
    if set(a)<={'1'}:return {k:a.get('1',0)*v for k,v in b.items()}
    if set(b)<={'1'}:return {k:b.get('1',0)*v for k,v in a.items()}
  raise ValueError(text)
 return repr(sorted(norm(tree).items()))

def argument(text,cpp=False):
 text=text.strip();strings=[]
 # Adjacent literals after Native // or C++ lexical concatenation have one value.
 while m:=re.search('('+literal+r')\s*(?://\s*)?('+literal+')',text):
  a,b=m[1],m[2];text=text[:m.start()]+json.dumps(a[1:-1]+b[1:-1])+text[m.end():]
 def stash(m):
  token='LITERAL_'+str(len(strings));strings.append((token,m[0][1:-1]));return token
 text=re.sub(literal,stash,text)
 text=text.replace('%','.').replace('.true.','true').replace('.false.','false').replace('Kokkos::','')
 # Existing source bool fields are uint8 in the mirror; the host print casts to bool.
 if cpp:
  text=re.sub(r'\bbool\(([^()]*)\)',r'\1',text)
  text=re.sub(r'\bdump_values\((.*)\)$',r'\1',text)
  while m:=re.search(r'\bsubview\(',text):
   raw,end=balanced(text,m.end());items=split_top(raw);base=items[0];indices=[]
   for value in items[1:]:
    if value=='ALL':indices.append(':')
    elif value.strip('() ') in param_indices:indices.append(value.strip('() '))
    else:indices.append(polynomial('('+value+')+1'))
   value=base if all(i==':' for i in indices) else base+'['+','.join(indices)+']'
   text=text[:m.start()]+value+text[end:]
 else:
  for m in reversed(list(re.finditer(r'\b([\w.]+)\(([^()]*)\)',text))):
   base=m[1];items=split_top(m[2])
   if not any(':' in item for item in items):continue
   indices=[':' if value==':' else value if value in param_indices else polynomial(value) for value in items]
   value=base if all(i==':' for i in indices) else base+'['+','.join(indices)+']'
   text=text[:m.start()]+value+text[m.end():]
 text=re.sub(r'[\s()]','',text).lower()
 for token,value in strings:text=text.replace(token.lower(),'STRING('+repr(value)+')')
 return text

calls_checked=0;arguments_checked=0;call_counts={}
callnames={'check_nan','check_negative','initialize_error_headers','check_parameters_api'}
for n in new:
 expected=[];actual=[]
 for st in statements(nb[n]):
  if m:=re.search(r'\bcall\s+(\w+)\s*(.*)',st,re.I):
   name=m[1];assert name in callnames,(n,name)
   raw=m[2].strip();aa=split_top(raw[1:-1]) if raw else []
   expected.append((name,[argument(a) for a in aa]))
 for m in re.finditer(r'(?<!\w)(?:f2k::)?('+'|'.join(callnames)+r')\s*\(',cb[n]):
  raw,end=balanced(cb[n],m.end());actual.append((m[1],[argument(a,True) for a in split_top(raw)] if raw.strip() else []))
 assert actual==expected,(n,'calls',[(a,e) for a,e in zip(actual,expected) if a!=e])
 call_counts[n]=len(expected);calls_checked+=len(expected);arguments_checked+=sum(len(aa) for name,aa in expected)

# Five stats-gated checks retain each exact field name and source order.
pdf=cb['pdf_closure_check'];native_pdf='\n'.join(statements(nb['pdf_closure_check']))
expected=re.findall(r'var_on_stats_list\(\s*stats,\s*"([^"]*)"',native_pdf)
actual=re.findall(r'var_on_stats_list\(\s*stats,\s*"([^"]*)"',pdf)
assert expected==actual==['wp4','wprtp2','wpthlp2','rcp2','wprtpthlp']
assert '"wpsclrtlp"' in pdf and '"pdf_params%w_1(1,:)"' in pdf

# Every newly ported WRITE payload, including values and literal concatenations.
# Existing four implementations are protected byte-for-byte above.
diagnostic_counts={};payload_count=0
for n in new:
 expected=[];actual=[]
 for st in statements(nb[n]):
  if m:=re.search(r'\bwrite\s*\(\s*fstderr\s*,\s*\*\s*\)\s*(.*)',st,re.I):
   if 'NaN in a hydrometeor model array' in m[1]:continue # Checked explicitly below.
   expected.append([argument(a) for a in split_top(m[1])])
 for m in re.finditer(r'diagnostic_message\s*<<\s*(.*?);',cb[n],re.S):
  if m[1].startswith('std::setprecision'):continue
  payload=[a.strip() for a in split_top(m[1],'<<') if a.strip() not in ['" "','"\\n"']]
  actual.append([argument(a,True) for a in payload])
 assert actual==expected,(n,'diagnostics',[(i,a,e) for i,(a,e) in enumerate(zip(actual,expected)) if a!=e],len(actual),len(expected))
 diagnostic_counts[n]=len(expected);payload_count+=sum(map(len,expected))
assert 'names_host(character, i - 1)' in cb['invalid_model_arrays']
assert 'character < 10' in cb['invalid_model_arrays']
assert '"NaN in a hydrometeor model array " + name' in cb['invalid_model_arrays']
assert "find_last_not_of(' ')" in cb['invalid_model_arrays']

# Compare host branch predicates independently after removing execution lambdas.
settings_native='\n'.join(statements(nb['check_clubb_settings_api']))
settings_cpp=cb['check_clubb_settings_api'];relations=[];native_conditions=[]
for st in statements(nb['check_clubb_settings_api']):
 if m:=re.match(r'(?:else\s*)?if\s*\(',st,re.I):
  raw,end=balanced(st,m.end());raw=raw.strip()
  if match:=re.search(r'\bany\s*\(',raw,re.I):
   inner,last=balanced(raw,match.end());relations.append(inner)
   raw=raw[:match.start()]+'(any_parameter_relation_'+str(len(relations))+' != 0)'+raw[last:]
  native_conditions.append(raw)
host_cpp=settings_cpp;kernel_conditions=[]
for m in reversed(list(re.finditer(r'KOKKOS_LAMBDA\([^{}]*?\)\s*\{',settings_cpp))):
 kernel,end=balanced(settings_cpp,m.end(),'{','}');start=m.start()
 condition=re.search(r'if\s*\(',kernel);assert condition
 raw,last=balanced(kernel,condition.end());kernel_conditions.insert(0,raw)
 host_cpp=host_cpp[:start]+' '*len(host_cpp[start:end])+host_cpp[end:]
cpp_conditions=[]
for m in re.finditer(r'\bif\s*\(',host_cpp):raw,end=balanced(host_cpp,m.end());cpp_conditions.append(raw)

def condition(text,cpp=False):
 if not cpp:
  for a,b in [('.and.','&&'),('.or.','||'),('.not.','!'),('/=','!=')]:text=re.sub(re.escape(a),b,text,flags=re.I)
 text=text.replace('f2k::','')
 # Shared metadata access replaces only the six source global sponge components.
 text=re.sub(r'environment\.(\w+_sponge_damp_settings)__l_sponge_damping',r'\1.l_sponge_damping',text)
 return argument(text,cpp)
assert [condition(a,True) for a in cpp_conditions]==[condition(e) for e in native_conditions],('settings conditions',[(a,e) for a,e in zip(cpp_conditions,native_conditions) if condition(a,True)!=condition(e)],len(cpp_conditions),len(native_conditions))
assert len(relations)==len(kernel_conditions)==11
# Reduction expressions use params(:,idx) pointwise at Native column+1.
for expected,actual in zip(relations,kernel_conditions):
 expected=expected.replace('%','.');expected=re.sub(r'params\(\s*:\s*,\s*(\w+)\s*\)',r'params(column,\1)',expected)
 actual=re.sub(r'params\(\s*\(1\s*-\s*1\s*\+\s*column\s*\*\s*\(1\)\)\s*,',r'params(column,',actual)
 assert condition(actual,True)==condition(expected),(expected,actual)

# Fatal updates/returns remain in exact host control flow. Native warnings whose
# fatal statements are commented out must stay warnings, including both-false.
for n in new:
 native='\n'.join(statements(nb[n]));cpp=cb[n]
 expected=len(re.findall(r'^err_info%err_code\s*=\s*clubb_fatal_error',native,re.M))
 actual=len(re.findall(r'deep_copy\(err_info.err_code,\s*clubb_fatal_error\)',cpp))
 assert actual==expected,(n,'fatal',expected,actual)
 expected=len(re.findall(r'^return\b',native,re.M));actual=len(re.findall(r'\breturn\b',cpp))
 assert actual==expected,(n,'return count',expected,actual)
assert settings_cpp.index('initialize_error_headers()')<settings_cpp.index('check_parameters_api(')<settings_cpp.index('driver_environment()')
assert '#ifdef CLUBB_GPU' in settings_cpp and '#endif' in settings_cpp
assert 'throw std::runtime_error("l_input_fields = .true. not usable when running on GPUs")' in settings_cpp
assert not re.search(r'\b(fma|__dmul_rn)\s*\(',target)

# Intrinsic/helper lowering: finite probes never stage valid profiles. Public
# UDTs and Native calls never cross a device kernel boundary.
cpp=clean(target);kernels=0
for m in re.finditer(r'KOKKOS_LAMBDA\([^{}]*?\)\s*\{',cpp):
 body,end=balanced(cpp,m.end(),'{','}');kernels+=1
 assert not re.search(r'\b(?:gr|err_info|pdf_params|clubb_config_flags|environment|sclr_idx)\.',body)
 assert 'f2k::' not in body and 'to_host_space(' not in body and 'std::cerr' not in body
 assert 'Kokkos::parallel_' not in body
# Ordered full-source nonfinite probes preserve hydrometeor/scalar loop order.
invalid_native='\n'.join(statements(nb['invalid_model_arrays']));invalid_cpp=cb['invalid_model_arrays']
expected=[];actual=[]
for m in re.finditer(r'any\(\s*\.not\.\s*ieee_is_finite\(',invalid_native,re.I):
 raw,end=balanced(invalid_native,m.end());expected.append(argument(raw))
for m in re.finditer(r'\bany_nonfinite\(',invalid_cpp):
 raw,end=balanced(invalid_cpp,m.end());actual.append(argument(raw,True))
assert actual==expected,(actual,expected)
nonfinite_probe_count=len(actual)
# Formula associations are exact for the two newly computed real quantities.
for n,lhs in [('rad_check','rvm'),('calculate_spurious_source','spurious_source')]:
 native=next(st.split('=',1)[1].strip() for st in statements(nb[n]) if re.match(lhs+r'\s*=',st,re.I))
 cpp=re.search(lhs+r'(?:\([^;]*?\))?\s*=\s*([^;]*);',cb[n])[1]
 if n=='rad_check':cpp=re.sub(r'\(slice_0,\s*slice_1\)','',cpp)
 assert ast.dump(ast.parse(native,mode='eval'))==ast.dump(ast.parse(cpp.strip(),mode='eval')),(n,native,cpp)

for filename in ['numerical_check.cpp','numerical_check.hpp']:
 text=(folder/filename).read_text();probe=mask(text)
 for left,right in [('(',')'),('[',']'),('{','}')]:
  depth=0
  for ch in probe:depth+=(ch==left)-(ch==right);assert depth>=0,(filename,left)
  assert depth==0,(filename,left,depth)
 assert all(line.rstrip()==line and '\t' not in line for line in text.splitlines()),filename
 # The four reviewed live blocks are deliberately not reformatted.
 cleaned=text
 if filename.endswith('.cpp'):
  for n in live:cleaned=cleaned.replace(definition(text,n)[2],'')
 assert all(len(line)<=100 for line in cleaned.splitlines()),filename
assert 'extern "C"' not in cpp and not re.search(r'\bf2k_proc_\w+',cpp)

# Match enclosing source IF branches for every new call, WRITE, fatal update and
# RETURN. This checks placement, beyond the independent linear payload scan.
def native_events(name):
 stack=[];result=[];relation=0
 for st in statements(nb[name]):
  inline=None
  if m:=re.match(r'(else\s*)?if\s*\(',st,re.I):
   cond,end=balanced(st,m.end());tail=st[end:].strip()
   if match:=re.search(r'\bany\s*\(',cond,re.I):
    if name=='check_clubb_settings_api':
     inner,last=balanced(cond,match.end());relation+=1
     cond=cond[:match.start()]+'(any_parameter_relation_'+str(relation)+' != 0)'+cond[last:]
   if m[1]:assert stack;stack.pop()
   cond=condition(cond)
   if tail.lower()=='then':stack.append(cond);continue
   inline=cond;st=tail
  elif re.match(r'^else\s*$',st,re.I):stack[-1]='ELSE('+stack[-1]+')';continue
  elif re.match(r'^end\s*if\b',st,re.I):assert stack;stack.pop();continue
  kind=None
  if re.search(r'^call\s+\w+',st,re.I):kind='CALL'
  elif re.search(r'^write\s*\(',st,re.I):kind='WRITE'
  elif re.search(r'^err_info%err_code\s*=\s*clubb_fatal_error',st,re.I):kind='FATAL'
  elif re.search(r'^return\b',st,re.I):kind='RETURN'
  elif re.search(r'^error\s+stop\b',st,re.I):kind='STOP'
  if kind:result.append((kind,tuple(stack+([inline] if inline else []))))
 assert not stack,name
 return result

def cpp_events(name):
 text=cb[name];probe=mask(text)
 for m in reversed(list(re.finditer(r'KOKKOS_LAMBDA\([^{}]*?\)\s*\{',probe))):
  body,end=balanced(text,m.end(),'{','}');text=text[:m.start()]+' '*(end-m.start())+text[end:]
 ranges=[]
 for m in re.finditer(r'\bif\s*\(',mask(text)):
  cond,end=balanced(text,m.end());start=end
  while text[start].isspace():start+=1
  assert text[start]=='{',(name,cond,text[start:start+60])
  body,end=balanced(text,start+1,'{','}');ranges.append((start,end,condition(cond,True)))
  after=end
  while after<len(text) and text[after].isspace():after+=1
  if text[after:after+4]=='else':
   after+=4
   while text[after].isspace():after+=1
   if text[after]=='{':
    body,last=balanced(text,after+1,'{','}')
    ranges.append((after,last,'ELSE('+condition(cond,True)+')'))
 patterns=[('CALL',r'(?<!\w)(?:f2k::)?(?:check_nan|check_negative|initialize_error_headers|check_parameters_api)\s*\('),
           ('WRITE',r'diagnostic_message\s*<<\s*(?!std::setprecision)'),
           ('FATAL',r'Kokkos::deep_copy\(err_info.err_code,\s*clubb_fatal_error\)'),
           ('RETURN',r'\breturn\b'),('STOP',r'\bthrow\s+std::runtime_error\(')]
 events=[]
 for kind,pattern in patterns:
  for m in re.finditer(pattern,text):
   # Negative lookahead can backtrack into whitespace; exclude precision directly.
   if kind=='WRITE' and text[m.end():].lstrip().startswith('std::setprecision'):continue
   events.append((m.start(),kind))
 if name=='invalid_model_arrays':
  events.append((text.index('std::cerr << "NaN in a hydrometeor model array'), 'WRITE'))
 result=[]
 for location,kind in sorted(events):
  guards=[guard for start,end,guard in sorted(ranges) if start<=location<end]
  result.append((kind,tuple(guards)))
 return result
# Intrinsics' syntax differs at this exact lowering boundary.
previous_condition=condition
def condition(text,cpp=False):
 if not cpp:
  text=re.sub(r'any\(\s*\.not\.\s*ieee_is_finite\(\s*(.*?)\s*\)\s*\)',r'any_nonfinite(\1)',text,flags=re.I)
  text=re.sub(r'ieee_is_finite\s*\(',r'isfinite(',text,flags=re.I)
 return previous_condition(text,cpp)
event_counts={}
for n in new:
 expected=native_events(n);actual=cpp_events(n)
 assert actual==expected,(n,'branch contexts',len(actual),len(expected),[(i,a,e) for i,(a,e) in enumerate(zip(actual,expected)) if a!=e][:8])
 event_counts[n]=len(expected)

record={
 'status':'Complete supported source draft; uncompiled and not runtime validated',
 'source':'src/CLUBB_core/numerical_check.F90; all1889lines read',
 'jax_scope':'All466lines read; full Native checks restored, not JAX reduced settings/schema',
 'routine_order':names,'argument_counts':counts,'new_procedure_count':len(new),
 'live_definition_byte_equality':preserved,'live_public_declaration_byte_equal':True,
 'new_ordered_call_counts':call_counts,'new_calls_checked':calls_checked,'new_arguments_checked':arguments_checked,
 'stats_gated_fields':expected if False else ['wp4','wprtp2','wpthlp2','rcp2','wprtpthlp'],
 'new_complete_diagnostic_payload_counts':diagnostic_counts,'diagnostic_payload_values':payload_count,
 'hydrometeor_name_width10_storage_and_trim_manually_checked':True,
 'settings_host_branch_predicates':len(native_conditions),'settings_device_ANY_relations':len(relations),
 'ordered_model_nonfinite_probes':nonfinite_probe_count,'device_kernel_bodies':kernels,
 'source_fatal_assignment_and_return_counts_match':True,'ordered_new_call_write_fatal_return_branch_events':event_counts,
 'manual_review':['All8new Native bodies reread end-to-end including diagnostics, complete comments and private generic dispatch', 'PDF first-column fields and typo wpsclrtlp retained; five exact stats gates', 'All source GPU warnings/ERRORSTOP and debug0 checks retained', 'Both-false relaxed-clipping source comments stay nonfatal; tau10 warnings keep >one OR <one', 'Eight host implementation checks accumulate general fatal flags without early return', 'NaN-before-negative/beginning-of-step/cold-reset behavior lives in unchanged definitions', 'Concrete source branches and all Native argument groups retained, no omitted option'],
 'adaptations':['Execution-space const LayoutStride profile Views preserve per-column strided callers', 'Required scalar surface API plus actual HostSpace vector adaptation through existing array owner', 'Private finite ANY reduction/dump formatter lower intrinsics and list-directed I/O only', 'calculate_spurious_source header placement enables native scalar device calls', 'Six source sponge globals read once through existing driver_environment metadata owner', 'GPU ERRORSTOP becomes host runtime_error, diagnostics use full-precision C++ streams'],
 'runtime':'No builds/tests/numerical execution; CPU BFB remains root integration gate',
 'sha256':{name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in ['numerical_check.cpp','numerical_check.hpp','f2k_dependencies.hpp']},
 'native_sha256':hashlib.sha256(source.encode()).hexdigest(),
 'jax_sha256':hashlib.sha256((root/'clubb_jax/src/CLUBB_core/numerical_check.py').read_bytes()).hexdigest()
}
(folder/'source_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ['status','argument_counts','new_calls_checked','new_arguments_checked','diagnostic_payload_values','settings_host_branch_predicates','settings_device_ANY_relations','device_kernel_bodies']},indent=2))
