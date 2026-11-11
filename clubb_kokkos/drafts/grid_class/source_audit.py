"""Read-only draft contract/formula/ownership audit; no model code is executed."""
from pathlib import Path
import hashlib,json,re,sys
D=Path(__file__).resolve().parent;R=D.parents[2]
sys.path.insert(0,str(R/'clubb_kokkos/src/f2k'))
from fortran_procedures import read_procedure,statements
from generate_derived_types import split_top
S=(R/'src/CLUBB_core/grid_class.F90').read_text();A=(R/'src/CLUBB_core/clubb_api_module.F90').read_text()
strip=re.compile(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"',re.S)
def clean(t):return strip.sub(lambda m:re.sub(r'[^\n]',' ',m[0]) if m[0].startswith('/') else m[0],t)
files=['prefix_definitions.cpp.inc','calc_zt2zm_weights.cpp.inc','calc_zm2zt_weights.cpp.inc','private_declarations.cpp.inc','procedure_declarations.hpp.inc','procedure_includes.hpp.inc','clubb_api_module.hpp.inc','api_includes.hpp.inc','cpp_includes.inc']
C={f:clean((D/f).read_text()) for f in files}
canonical=['setup_grid','cleanup_grid_api','setup_grid_heights','calc_zt2zm_weights','calc_zm2zt_weights'];defs=[]
for f in files[:3]:
 for m in re.finditer(r'\b(?:static\s+)?void\s+(\w+)\s*\((.*?)\)\s*\{',C[f],re.S):
  n=m[1];args=[x.strip().split()[-1] for x in split_top(m[2])]
  assert args==[a.name for a in read_procedure(R/'src/CLUBB_core/grid_class.F90',n).arguments],(n,args)
  defs.append(n)
assert defs==canonical
B=(D/files[0]).read_text();b=C[files[0]]
# Every allocated/released component retains the actual Native order.
alloc=['zm','zt','dzm','dzt','invrs_dzm','invrs_dzt','weights_zm2zt','weights_zt2zm']
assert re.findall(r'if\s*\(gr\.(\w+)\.is_allocated\(\)\)',b)==alloc
assert re.findall(r'gr\.(\w+)\s*=\s*\{\};',b)==alloc
assert b.count('Kokkos::WithoutInitializing')==8
assert b.index('calc_zm2zt_weights(nzt')<b.index('calc_zt2zm_weights(nzm')<b.index('gr.k_lb_zm = 1;')
assert 'sfc_elevation(i) > zm(i, 0)' in b and 'Kokkos::Min<int>(first_bad_column)' in b
assert 'Kokkos::subview(err_info.err_code, first_bad_column)' in b
assert b.count('Kokkos::deep_copy(err_info.err_code, clubb_fatal_error)')==6
assert b.index('if (!l_implemented && grid_type != 1')<b.index('const auto zm = gr.zm;',b.index('void setup_grid_heights'))
# Eight exact source weight expressions after direct one-based-to-zero-based translation.
w1=re.sub(r'\s+','',C[files[1]]);w2=re.sub(r'\s+','',C[files[2]])
expected=[
'weights_zt2zm(i,k,t_above)=(zm(i,k)-zt(i,k-1))/(zt(i,k)-zt(i,k-1));',
'weights_zt2zm(i,k,t_below)=(zt(i,k)-zm(i,k))/(zt(i,k)-zt(i,k-1));',
'weights_zt2zm(i,0,t_above)=(zm(i,0)-zt(i,0))/(zt(i,1)-zt(i,0));',
'weights_zt2zm(i,0,t_below)=(zt(i,1)-zm(i,0))/(zt(i,1)-zt(i,0));',
'weights_zt2zm(i,nzm-1,t_above)=(zm(i,nzm-1)-zt(i,nzt-2))/(zt(i,nzt-1)-zt(i,nzt-2));',
'weights_zt2zm(i,nzm-1,t_below)=(zt(i,nzt-1)-zm(i,nzm-1))/(zt(i,nzt-1)-zt(i,nzt-2));',
'weights_zm2zt(i,k,m_above)=(zt(i,k)-zm(i,k))/(zm(i,k+1)-zm(i,k));',
'weights_zm2zt(i,k,m_below)=(zm(i,k+1)-zt(i,k))/(zm(i,k+1)-zm(i,k));']
for e in expected:assert e in (w1 if 'zt2zm' in e else w2),e
# All four public API specifics keep exact names/order, exposed by Native generic overloads.
api_names=['setup_grid_heights_api_single_col','setup_grid_heights_api_multi_col','setup_grid_api_single_col','setup_grid_api_multi_col']
matches=list(re.finditer(r'void\s+(setup_grid(?:_heights)?_api)\s*\((.*?)\)\s*\{',C['clubb_api_module.hpp.inc'],re.S));assert len(matches)==4
for n,m in zip(api_names,matches):assert [x.strip().split()[-1] for x in split_top(m[2])]==[a.name for a in read_procedure(R/'src/CLUBB_core/clubb_api_module.F90',n).arguments],n
api=(D/'clubb_api_module.hpp.inc').read_text();assert api.count('Error in CLUBB calling setup_grid_heights\\n')==2 and api.count('Error in CLUBB calling setup_grid\\n')==2
# No owning UDT is captured; view components/metadata are bound outside kernels.
kernels=[]
for f in files[:3]:
 t=C[f]
 for m in re.finditer(r'KOKKOS_LAMBDA\s*\([^)]*\)\s*\{',t):
  i=m.end();depth=1
  while depth:depth+=(t[i]=='{')-(t[i]=='}');i+=1
  body=t[m.end():i-1];assert not re.search(r'\b(?:gr|err_info)\s*\.|f2k::|std::cerr|\bfma\s*\(',body),(f,body);kernels.append(f)
assert len(kernels)==5
# Native descriptions/weight derivations are retained verbatim as comment text.
for n,f in [('calc_zt2zm_weights',files[1]),('calc_zm2zt_weights',files[2])]:
 native=re.search(r'(?ims)^  subroutine '+n+r'\b.*?^  end subroutine '+n+r'\b',S)[0];comment=re.search(r'(?ms)! Description:.*?^\s*!\s*-{20,}\s*$',native)[0]
 for l in comment.splitlines():assert l.split('!',1)[-1].strip() in (D/f).read_text(),(n,l)
# Separate cleanup audit checks exact executable tokens before/after formatting.
lexer=re.compile(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|\w+|[^\s]')
def tokens(s):return lexer.findall(re.sub(r'//[^\n]*|/\*.*?\*/','',s,flags=re.S))
cleanup=[]
for f in files:
 p=D/f;before=Path('/tmp/grid_before_'+f)
 if before.exists():assert tokens(before.read_text())==tokens(p.read_text());cleanup.append({'file':f,'exact_executable_token_equality':True,'tokens':len(tokens(p.read_text()))})
 assert all(not l.rstrip('\n').endswith((' ','\t')) for l in p.read_text().splitlines(True)),f
 long=[i for i,l in enumerate(p.read_text().splitlines(),1) if len(l)>100];assert not long,(f,long)
hashes={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in files}
source={'review_kind':'Owner full five-body/source-contract audit, separate cleanup. Not an independent review.','Native_source_sha256':hashlib.sha256(S.encode()).hexdigest(),'hashes':hashes,'new_definitions_in_Native_merge_order':defs,'argument_signatures':5,'Native_API_specifics':api_names,'API_generic_overloads':4,'weight_formulas':8,'View_only_kernels':len(kernels),'allocation_and_cleanup_order':alloc,'reviewed_branches':['implemented/standalone','uniform/stretched_zt/stretched_zm','ascending/descending','all four first-column range failure exits','allocation failure and already-allocated components','debug-gated height fatal handling','first failing surface column only','direct setup_heights invalid grid type before writes','signed spacings and direct reciprocals','two boundary linear extensions','direction metadata only after weights'],'deliberate_language_boundaries':['read_grid_heights remains entire existing Native file/validation service','Kokkos reference-count release cannot return Fortran DEALLOCATE STAT','bad-allocation exception handling does not assert recovery from a Kokkos runtime abort','API single-column packing and generic array staging use the shared owner','independent boundary/interior weights share one kernel; no arithmetic reassociation or FMA added'],'runtime_status':'Draft only; no shared edits/build/tests/CPU/GPU/BFB claim. All runtime branches pending root integration.'}
(D/'SOURCE_AUDIT.json').write_text(json.dumps(source,indent=2)+'\n');(D/'CLEANUP_AUDIT.json').write_text(json.dumps({'hashes':hashes,'checks':cleanup,'max_columns':100,'no_trailing_whitespace':True,'source_grouping':'Native physical argument/call groups retained, splitting typed C++ width only.','runtime_status':source['runtime_status']},indent=2)+'\n');print('Grid source/cleanup audits passed: five bodies, four API overloads, eight weight formulas, five kernels.')
