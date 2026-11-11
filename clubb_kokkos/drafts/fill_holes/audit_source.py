"""Separate source/interface/ownership audit; no compilation or runtime model test."""
from pathlib import Path
import hashlib,json,re
folder=Path(__file__).resolve().parent;root=folder.parents[2]
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

outline=json.loads((folder/'SOURCE_OUTLINE.json').read_text())
supported=[r for r in outline['routines'] if r['supported']]
source=(root/'src/CLUBB_core/fill_holes.F90').read_text()
cpp=(folder/'fill_holes.cpp').read_text();hpp=(folder/'fill_holes.hpp').read_text()
def strip(s):
    return re.sub(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'',lambda m:re.sub(r'[^\n]',' ',m[0]) if m[0].startswith(('/',)) else m[0],s,flags=re.S)
c=strip(cpp);h=strip(hpp);definitions={};order=[]
for m in re.finditer(r'\bvoid\s+(\w+)\s*\(',c):
    arg,end=balanced(c,m.end()-1)
    if not c[end:].lstrip().startswith('{'):continue
    name=m[1]
    if name not in [r['name'] for r in supported]:continue
    arg=re.sub(r'Kokkos::View<[^>]+>','Char2D',arg)
    args=[re.search(r'(\w+)\s*$',a)[1] for a in split(arg)]
    brace=c.index('{',end);body,_=balanced(c,brace,'{','}')
    definitions[name]={'arguments':args,'body':body};order.append(name)
assert order==[r['name'] for r in supported],order
for r in supported:
    expected=[a['name'] for a in r['arguments']]
    assert definitions[r['name']]['arguments']==expected,r['name']
    if r['name'] in ('fill_holes_global','fill_holes_sliding_window'):continue
    m=re.search(r'\bvoid\s+'+r['name']+r'\s*\(',h);arg,end=balanced(h,m.end()-1)
    arg=re.sub(r'Kokkos::View<[^>]+>','Char2D',arg)
    actual=[re.search(r'(\w+)\s*(?:=\s*nullptr)?\s*$',a)[1] for a in split(arg)]
    assert actual==expected,(r['name'],actual,expected)
    assert h[end:].lstrip().startswith(';'),r['name']
    body=cpp[cpp.index('void '+r['name']):]
    assert '// Description:' in body[:body.find('\n}',body.find('{'))]
device_f2k=[];owning_capture=[];kernels=0
for src in (c,h):
    for m in re.finditer(r'KOKKOS_LAMBDA\s*\(',src):
        args,end=balanced(src,m.end()-1);brace=src.index('{',end);body,_=balanced(src,brace,'{','}');kernels+=1
        if 'f2k::' in body:device_f2k.append(body)
        if re.search(r'\b(?:gr|hm_metadata|stats)\s*\.',body):owning_capture.append(body)
assert not device_f2k and not owning_capture
assert len(re.findall(r'std::string\s*\*\s*name_\w+\s*=\s*nullptr',hpp))==5
assert 'bool l_mfl' not in cpp
assert 'std::exit(EXIT_SUCCESS)' in cpp and 'std::exit(EXIT_FAILURE)' in cpp
assert 'grid_dir_indx * (k - upper_hf_level)' in cpp
assert re.search(r'Kokkos::parallel_for\(\s*"fill_holes_sliding_window::columns"',cpp)
assert 'widening_windows(' not in cpp and 'smart_window(' not in cpp and 'fill_holes_parallel(' not in cpp
assert 'if (stats.flags().l_sample)' in cpp
hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'src/CLUBB_core/fill_holes.F90',root/'clubb_jax/src/CLUBB_core/fill_holes.py',folder/'fill_holes.cpp',folder/'fill_holes.hpp',folder/'f2k_dependencies.hpp']}
manual=[
 'Global: source-order numerator/denominator/clipped numerator, original any-hole gate and eps/2 guard, ordered MIN/MAX operands, normalization and descending-range bounds reread.',
 'Sliding: all-field entry/fallback ANY, original start/stop and denominator j bounds, Native direction-dependent overlapping k recurrence, exact average->clip->average->scale sequence reread.',
 'TKE: early ANY return, exact Native sufficient/insufficient/no-donor predicates and assignments, source epsilon kind and debug7warning conditions/payloads reread. Device science emitted from Native; host diagnostics use pre-update snapshots.',
 'One-level frozen fill: strict negative-vs-else totals in species order, epsilon no-material return, insufficient/complete formulas, warning multiplicity and conservation/negative-result checks reread.',
 'Frozen parent: same exact compact species subset in Native order, call only when a level has negative frozen mixing ratios, correct packed/unpacked output and unchanged nonfrozen fields.',
 'Vapor: partial vapor mutation before unknown-species ERROR STOP retained, unknown/no-hole path returns; strict hole threshold, rain/frozen heat expressions and assignment order reread.',
 'Driver:8 Native ordered budget sites mapped through local name/species/phase adapter; sampling/any-hole/mix-ratio flags, species and column loops, conditional vertical/vapor calls, small-mass vapor/heat compensation, final concentration copy reread.',
 'Concentration: per-species public metadata functions at host boundary, Native1based indices, ((radius*radius)*radius), distinct liquid/ice densities, strict positive mixing-ratio branch and unchanged mass fields reread.',
 'Stats names:9 exact named species cases, fixed-width trailing-blank semantics, optional5outputs/default absence, reset speeds/blank names and same velocities reread.',
 'Ownership/format:8public+2private definitions in Native order, interfaces/intents/optionals, ordinary groups, Native descriptions, field/axis units, per-call storage, read/write loops, source-only math adaptations, and complete unsupported scope reread.'
]
result={'status':'Source audit/cleanup complete for isolated draft. Uncompiled and numerically unverified.',
 'supported_definitions':len(supported),'public_procedures':8,'excluded_native_private_methods':[r['name'] for r in outline['routines'] if not r['supported']],
 'signature_name_orders_match':True,'public_header_signatures_and_prototype_terminators_match':True,
 'optional_stats_name_outputs':5,'device_lambda_bodies':kernels,'device_f2k_calls':device_f2k,'owning_udt_captures':owning_capture,
 'manual_source_reviews':manual,'hashes':hashes,'builds_tests_or_numerical_claims':False,
 'remaining':'Root integration/public lower bridge generation/constants and CPU whole-model BFB; performance of per-column microphysics API packing is unmeasured.'}
(folder/'source_audit.json').write_text(json.dumps(result,indent=2)+'\n')
cleanup={'separate_source_cleanup':True,'findings_corrected':['Prototype semicolon originally followed intent // and was swallowed; now standalone syntax.','Quoted diagnostic commas used declaration splitter; fixed quote-aware rendering preserving exact labels.','Explicit device-template selection prevents host wrapper recursion on compatible rank1 Serial Views.','Public vertical/TKE generic adapters now use actual existing array owner for strided caller sections.','Native setup_stats_names species-case rationale retained; supplied missing argument/unit descriptions.'],
 'cpp_max_line_width':max(map(len,cpp.splitlines())),'hpp_max_line_width':max(map(len,hpp.splitlines())),
 'source_or_runtime_tests':'Source inspection only; no build, model run, unit test, numerical or performance result.'}
(folder/'cleanup_audit.json').write_text(json.dumps(cleanup,indent=2)+'\n')
print(json.dumps({'definitions':len(supported),'kernels':kernels,'hashes':hashes},indent=2))
