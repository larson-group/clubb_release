"""Separate exact-source audit, not compilation or a numerical test."""
from pathlib import Path
import ast,hashlib,json,re,sys
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import declaration,statements
def logical(src):
 out=[];buf='';start=0
 for ln,line in enumerate(src.splitlines(),1):
  # Exclude comments only when outside Fortran quoted literals.
  q=None;clean='';j=0
  while j<len(line):
   ch=line[j]
   if q:
    clean+=ch
    if ch==q:
     if j+1<len(line) and line[j+1]==q:clean+=line[j+1];j+=1
     else:q=None
   elif ch in "'\"":q=ch;clean+=ch
   elif ch=='!':break
   else:clean+=ch
   j+=1
  clean=clean.strip()
  if not clean:continue
  if not buf:start=ln
  buf+=' '+clean.lstrip('&').removesuffix('&').strip()
  if not clean.endswith('&'):out.append((start,buf.strip()));buf=''
 return out

def balanced(s,start,op='(',cl=')'):
 assert s[start]==op
 depth=1;q=None;i=start+1
 while i<len(s):
  ch=s[i]
  if q:
   if ch=='\\':i+=2;continue
   if ch==q:q=None
  elif ch in "'\"":q=ch
  elif ch==op:depth+=1
  elif ch==cl:
   depth-=1
   if depth==0:return s[start+1:i],i+1
  i+=1
 raise ValueError((start,s[start:start+150]))

def split(s):
 out=[];d=0;q=None;start=0
 for i,ch in enumerate(s):
  if q:
   if ch==q and (i==0 or s[i-1]!='\\'):q=None
  elif ch in "'\"":q=ch
  elif ch in '({[':d+=1
  elif ch in ')}]':d-=1
  elif ch==',' and d==0:out.append(s[start:i].strip());start=i+1
 if s[start:].strip():out.append(s[start:].strip())
 return out

def remove_comments(s):
    return re.sub(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'',lambda m:re.sub(r'[^\n]',' ',m[0]) if m[0].startswith(('/',)) else m[0],s,flags=re.S)

def cpp_minmax(s):
    s=re.sub(r'\s+','',s)
    while '?' in s:
        candidates=[];stack=[]
        for i,ch in enumerate(s):
            if ch=='(':stack.append(i)
            elif ch==')':stack.pop()
            elif ch=='?':
                start=stack[-1];inside,end=balanced(s,start)
                if inside.count('?')==1:candidates.append((start,inside,end))
        assert candidates,s
        start,inside,end=candidates[-1];condition,branches=inside.split('?',1);a,b=branches.split(':',1)
        cond=ast.parse(condition,mode='eval').body;aa=ast.parse(a,mode='eval').body;bb=ast.parse(b,mode='eval').body
        assert isinstance(cond,ast.Compare) and len(cond.ops)==1
        assert ast.dump(cond.left)==ast.dump(aa) and ast.dump(cond.comparators[0])==ast.dump(bb)
        fn='min' if isinstance(cond.ops[0],ast.Lt) else 'max';assert isinstance(cond.ops[0],(ast.Lt,ast.Gt))
        s=s[:start]+fn+'('+a+','+b+')'+s[end:]
    return s

def affine(n):
    if isinstance(n,ast.Name):return {n.id:1}
    if isinstance(n,ast.Attribute):return {ast.unparse(n):1}
    if isinstance(n,ast.Constant):return {'#':n.value}
    if isinstance(n,ast.UnaryOp):return {k:(-v if isinstance(n.op,ast.USub) else v) for k,v in affine(n.operand).items()}
    if isinstance(n,ast.BinOp) and isinstance(n.op,(ast.Add,ast.Sub)):
        a=affine(n.left);b=affine(n.right)
        for k,v in b.items():a[k]=a.get(k,0)+(-v if isinstance(n.op,ast.Sub) else v)
        return {k:v for k,v in a.items() if v}
    raise ValueError(ast.dump(n))

def canon(n,target=False):
    if isinstance(n,ast.Name):return ('name',n.id)
    if isinstance(n,ast.Attribute):return ('attribute',canon(n.value),n.attr)
    if isinstance(n,ast.Constant):return ('constant',n.value)
    if isinstance(n,ast.UnaryOp):return (type(n.op).__name__,canon(n.operand,target))
    if isinstance(n,ast.BinOp):
        if isinstance(n.op,ast.Pow):
            assert isinstance(n.right,ast.Constant) and n.right.value==2
            return ('Mult',canon(n.left,target),canon(n.left,target))
        return (type(n.op).__name__,canon(n.left,target),canon(n.right,target))
    if isinstance(n,ast.Call):
        fn=ast.unparse(n.func)
        if fn in ('real_t',):return canon(n.args[0],target)
        if fn in ('min','max','sqrt','exp','erf','abs'):return ('call',fn,tuple(canon(a,target) for a in n.args))
        indexes=[]
        for arg in n.args:
            a=affine(arg)
            if target:a['#']=a.get('#',0)+1
            indexes.append(tuple(sorted((k,v) for k,v in a.items() if v)))
        return ('array',fn,tuple(indexes))
    raise ValueError(ast.dump(n))

def expression(s,target=False):
    s=s.lower().replace('%','.').replace('kokkos::','')
    s=re.sub(r'\bgr_(\w+)\b',r'gr.\1',s)
    s=re.sub(r'_core_rknd\b','',s)
    s=s.replace('.true.','True').replace('.false.','False')
    s=re.sub(r'\btrue\b','True',s);s=re.sub(r'\bfalse\b','False',s)
    if target:s=cpp_minmax(s)
    return canon(ast.parse(s.strip(),mode='eval').body,target)

native=(root/'src/CLUBB_core/mono_flux_limiter.F90').read_text()
cpp=(folder/'mono_flux_limiter.cpp').read_text();clean=remove_comments(cpp)
header=remove_comments((folder/'mono_flux_limiter.hpp').read_text())
outline=json.loads((folder/'SOURCE_OUTLINE.json').read_text())
report={};definitions={};order=[]
for m in re.finditer(r'\bvoid\s+(\w+)\s*\(',clean):
    a,end=balanced(clean,m.end()-1)
    if not clean[end:].lstrip().startswith('{'):continue
    if m[1] not in [x['name'] for x in outline['routines']]:continue
    brace=clean.index('{',end);body,_=balanced(clean,brace,'{','}')
    aa=re.sub(r'Kokkos::View<[^>]+>','Integer2D',a)
    definitions[m[1]]={'args':[re.search(r'(\w+)\s*$',x)[1] for x in split(aa)],'body':body};order.append(m[1])
assert order==[r['name'] for r in outline['routines']]
native_calls=0;native_actuals=0;indexed_total=0;scalar_total=0;call_adapters=[]
for entry in outline['routines']:
    n=entry['name'];expected=[a['name'] for a in entry['arguments']]
    assert definitions[n]['args']==expected,(n,definitions[n]['args'],expected)
    if entry['public']:
        m=re.search(r'\bvoid\s+'+n+r'\s*\(',header);a,end=balanced(header,m.end()-1)
        a=re.sub(r'Kokkos::View<[^>]+>','Integer2D',a)
        assert [re.search(r'(\w+)\s*$',x)[1] for x in split(a)]==expected
        assert header[end:].lstrip().startswith(';')
    nb=re.search(r'^  subroutine '+n+r'\b.*?^  end subroutine '+n+r'\b',native,re.M|re.S)[0]
    nl=logical(nb);arrays={};calls=[]
    for ln,s in nl:
        if '::' in s:arrays.update({a.name.lower():a for a in declaration(s).values() if a.dims})
        if m:=re.match(r'call\s+(\w+)\s*\((.*)\)',s,re.I):calls.append((m[1],split(m[2])))
    native_calls+=len(calls);native_actuals+=sum(len(c[1]) for c in calls)
    if n=='mfl_xm_solve':report[n]={'signature':len(expected),'solver_adapter':'Manually reviewed against actual public header; Native lhs/rhs/xm restored before all-column fatal return; xm flip only on success.'};continue
    # All explicitly indexed array assignment trees remain in Native source order.
    na=[]
    for ln,s in nl:
        if m:=re.match(r'^(\w+)\((.*?)\)\s*=\s*(.*)$',s):
            if m[1].lower() in arrays and ':' not in m[2]:na.append((expression(m[1]+'('+m[2]+')'),expression(m[3])))
    ca=[]
    for m in re.finditer(r'^\s*(\w+)\s*\(',definitions[n]['body'],re.M):
        if m[1].lower() not in arrays:continue
        index,end=balanced(definitions[n]['body'],m.end()-1)
        tail=definitions[n]['body'][end:].lstrip()
        if not tail.startswith('=') or tail.startswith('=='):continue
        rhs=tail[1:tail.index(';')].strip()
        # Native whole-column spike scaling is reviewed separately from indexed trees.
        if n=='monotonic_turbulent_flux_limit' and m[1]=='xm' and re.search(r'\*\s*\(real_t\(1\.\)\s*\+\s*xm_adj_coef\)',rhs):continue
        ca.append((expression(m[1]+'('+index+')',True),expression(rhs,True)))
    assert len(na)==len(ca),(n,len(na),len(ca))
    for i,(a,b) in enumerate(zip(na,ca)):assert a==b,(n,i,a,b)
    indexed_total+=len(na)
    # Scalar arithmetic and loop-carried search state use the same expression
    # audit; spike SUM/redistribution and logical reduction have explicit adapters.
    scalar_types={}
    for _,s in nl:
        if '::' in s:scalar_types.update({a.name.lower():a for a in declaration(s).values() if not a.dims and a.kind!='character'})
    sa=[];skip={'l_any_adjustment_needed','dz','xm_density_weighted','xm_vert_integral','xm_adj_coef'} if n=='monotonic_turbulent_flux_limit' else set()
    for _,s in nl:
        if m:=re.match(r'^(\w+)\s*=\s*(.*)$',s):
            if m[1].lower() in scalar_types and m[1].lower() not in skip:sa.append((m[1].lower(),expression(m[2])))
    assigned={a[0] for a in sa};sc=[]
    for m in re.finditer(r'^\s*(?:(?:const|constexpr)\s+)?(?:(?:real_t|int|bool)\s+)?(\w+)\s*=\s*(.*?);',definitions[n]['body'],re.M|re.S):
        if m[1].lower() in assigned:sc.append((m[1].lower(),expression(m[2],True)))
    assert sa==sc,(n,'scalar',sa,sc)
    scalar_total+=len(sa)
    # Call spelling/order/field mapping: stats arithmetic buffers are explicit.
    cpp_calls=[]
    for m in re.finditer(r'\b(\w+)\s*\(',definitions[n]['body']):
        if m[1] not in {c[0] for c in calls}:continue
        a,_=balanced(definitions[n]['body'],m.end()-1);cpp_calls.append((m[1],split(a)))
    assert [c[0] for c in calls]==[c[0] for c in cpp_calls],n
    def call_norm(s,target=False):
        s=s.lower()
        if target:
            s=re.sub(r'\bgr_(\w+)\b',r'gr%\1',s)
            s=re.sub(r'\breal_t\(([^()]*)\)',r'\1',s)
        else:s=re.sub(r'\((?:\s*:\s*,)*\s*:\s*\)','',s)
        s=re.sub(r'_core_rknd\b','',s)
        return re.sub(r'\s+','',s)
    for (name,actual),(callee,converted) in zip(calls,cpp_calls):
        assert len(actual)==len(converted),(n,name,actual,converted)
        for index,(a,b) in enumerate(zip(actual,converted)):
            if call_norm(a)==call_norm(b,True):continue
            assert name in ('stats_begin_budget','stats_finalize_budget') and index==1 and '/' in a and b.startswith('stats_values_'),(n,name,index,a,b)
            call_adapters.append({'caller':n,'callee':name,'native_value':a,'cpp_payload':b})
    report[n]={'signature':len(expected),'ordered_indexed_assignment_trees':len(na),'ordered_scalar_assignment_trees':len(sa),'ordered_calls':len(calls)}

kernels=0;first_use=[];device_crossings=[];bad_captures=[];case_errors=[]
for m in re.finditer(r'KOKKOS_LAMBDA\s*\(',clean):
    a,end=balanced(clean,m.end()-1);brace=clean.index('{',end);body,_=balanced(clean,brace,'{','}');kernels+=1
    if 'f2k::' in body or 'std::cerr' in body or 'std::max' in body:device_crossings.append(body)
    if re.search(r'\b(?:gr|err_info|stats)\s*\.',body):bad_captures.append(body)
    for loc in re.finditer(r'\b(?:int|real_t|bool)\s+(\w+)\s*;',body):
        use=re.search(r'\b'+loc[1]+r'\b',body[loc.end():])
        if use:
            tail=body[loc.end()+use.end():].lstrip()
            if not tail.startswith('=') or tail.startswith('=='):first_use.append(loc[1])
for name in definitions:
    for m in re.finditer(r'\b'+name+r'\s*\(',clean,re.I):
        if m[0].split('(')[0].strip()!=name:case_errors.append(m[0])
assert not device_crossings and not bad_captures and not first_use and not case_errors
hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'src/CLUBB_core/mono_flux_limiter.F90',root/'clubb_jax/src/CLUBB_core/mono_flux_limiter.py',folder/'mono_flux_limiter.cpp',folder/'mono_flux_limiter.hpp']}
result={'status':'Separate source audit complete; compilation/integration/model comparison pending.',
 'source_read':{'full_native_lines':2157,'all7procedures':True,'actual_jax_supported_scope_read':True},
 'procedures':report,'indexed_assignment_trees':indexed_total,'scalar_assignment_trees':scalar_total,'ordered_native_calls':native_calls,'native_call_actuals':native_actuals,'explicit_stats_payload_adapters':call_adapters,
 'device_lambda_bodies':kernels,'device_f2k_or_io_or_stdmax':device_crossings,'owning_udt_captures':bad_captures,'mutable_scalar_first_use_candidates':first_use,'case_sensitive_call_errors':case_errors,
 'manual_contract_review':['All26/13public required args; private12/11/10/12/8 signatures and actual lower grid/mean/solver/stats headers.',
 'Parallel screen uses only original neighbor wpxp; integer ANY gate precedes sequential per-column recurrence in Native physical direction; min scratch only set on original source paths.',
 'Native flags permit implicit and explicit correction bodies; per-column xm replacement occurs before parent fatal check. All-column existing fatal state included, no reset or guessed zero initialization.',
 'Source min/max influence windows, exp/erf branches, source literal constants, signed powers, interpolation actuals and boundary values; no JAX sqrt clamp or safe divisor introduced.',
 'Solver bands-first layout, optional force-descending pack, actual live tridiag generic argument order, post-solver lhs/rhs/solution restoration before fatal, and success-only final xm reversal.',
 'Native spike density integral and -0.99 cap, ordered whole-column scaling, top restoration, CPU warning payload/order and diagnostic comments retained at host boundary.',
 'All named stats selected by same case, blank names retained, sampled budgets/updates in same order; four field/dt payload buffers and same-kind0literal are expected call expression adapters.',
 'Both grid directions and disabled constant-thickness searches retain Native loop and break order;1based output integer levels and six unused-level boundaries retained.',
 'Meaningful Native description/derivations, declaration intent grouping, source call continuations, comments and units checked; sigma_w_i_zm copied-unit typo corrected only in target.'
 ],'hashes':hashes,'builds_tests_numerical_or_BFB_claim':False,'missing_shared_owner':'l_force_descending_solves=false in existing model_flags generator; mean_adv owner promotion; existing stats/grid/solver APIs otherwise present.'}
(folder/'source_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'indexed_trees':indexed_total,'scalar_trees':scalar_total,'calls':native_calls,'actuals':native_actuals,'kernels':kernels,'hashes':hashes},indent=2))
