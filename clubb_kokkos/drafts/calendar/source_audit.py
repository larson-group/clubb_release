"""Audit source contracts/expressions only; do not execute calendar calculations."""
from pathlib import Path
import re,sys,json,hashlib,ast
D=Path(__file__).resolve().parent;R=D.parents[2];P=R/'src/CLUBB_core/calendar.F90';S=P.read_text();T=(D/'calendar.cpp').read_text();C=re.sub(r'//[^\n]*','',T);sys.path.insert(0,str(R/'clubb_kokkos/src/f2k'));from fortran_procedures import read_procedure,statements
from generate_derived_types import split_top
helper=(R/'clubb_kokkos/drafts/advance_helper_module/source_audit.py').read_text();exec(helper[helper.index('class Cpp:'):helper.index('reviews=[]')]);Cpp.ops['%']=6
names=json.loads((D/'SOURCE_OUTLINE.json').read_text())['native_inventory_and_order'];definitions=[]
for m in re.finditer(r'\b(?:int|bool|void)\s+('+ '|'.join(names)+r')\s*\((.*?)\)\s*\{',C,re.S):
 e=m.end();depth=1
 while depth:depth+=(C[e]=='{')-(C[e]=='}');e+=1
 args=[a.strip().split()[-1] for a in split_top(m[2])];assert args==[a.name for a in read_procedure(P,m[1]).arguments];definitions.append((m[1],C[m.end():e-1]))
assert [n for n,b in definitions]==names

def normalize(x):
 if isinstance(x,tuple):
  if x[:2]==('call','mod'):return ('bin','%',*(normalize(a) for a in x[2]))
  if x[:2] in [('call','real'),('call','as_double')]:return normalize(x[2][0])
  if x[:2]==('call','int') and normalize(x[2][0])[:2]==('call','floor'):return normalize(x[2][0])
  if x[:2]==('call','std::floor'):return ('call','floor',normalize(x[2]))
  return tuple(normalize(y) for y in x)
 return x
CppOrigInit=Cpp.__init__
def init(self,t,arrays):return CppOrigInit(self,t.replace('static_cast<double>','as_double'),arrays)
Cpp.__init__=init
reviews=[]
for n,b in definitions:
 raw=re.search(r'^[ \t]*(?:(?:integer|logical)\s+function|subroutine)\s+'+n+r'\b.*?^[ \t]*end (?:function|subroutine)\s+'+n+r'\b',S,re.M|re.S|re.I).group()
 ff=[]
 for line in statements(raw):
  if n!='gregorian2julian_day' and (m:=re.match(r'([\w]+)\s*=\s*(.*)',line)):ff.append((m[1].lower(),normalize(fexpr(m[2]))))
 if n=='gregorian2julian_day':
  assert 'for (int m = 1; m <= month - 1; ++m)' in b and 'previous_month_days += days_per_month[m - 1];' in b
  assert 'int gregorian2julian_day = day + previous_month_days;' in b
  assert 'gregorian2julian_day = gregorian2julian_day + 1;' in b
  assert 'throw std::runtime_error("Problem with Julian day conversion in gregorian2julian_day.");' in b
  predicates=[normalize(fexpr(m[1])) for line in statements(raw) if (m:=re.match(r'if\s*\((.*)\)\s*then$',line,re.I))]
  cpp_predicates=[normalize(Cpp(m[1],set()).parse()) for m in re.finditer(r'if\s*\((.*?)\)\s*{',b,re.S)]
  assert predicates==cpp_predicates,(predicates,cpp_predicates)
  reviews.append({'name':n,'Native_ordered_month_SUM_and_function_result':True,'leap_correction_and_exact_fatal_literal':True});continue
 cc=[]
 for line in b.splitlines():
  if m:=re.match(r'\s*([\w]+)\s*=\s*(.*?);\s*$',line):cc.append((m[1].lower(),normalize(Cpp(m[2],set()).parse())))
 # Multiline assignment/return expressions retain their exact operator trees.
 if n in ['gregorian2julian_date','leap_year']:
  expr=re.search(r'\breturn\s+(.*?);',b,re.S)[1];cc.append((n,normalize(Cpp(expr,set()).parse())))
 if n=='compute_current_date_api':
  cc=[]
  for m in re.finditer(r'\b(days_since_1jan4713bc|days_since_start|seconds_since_current_date)\s*=\s*(.*?);',b,re.S):cc.append((m[1],normalize(Cpp(m[2],set()).parse())))
 assert ff==cc,(n,ff,cc)
 reviews.append({'name':n,'ordered_assignment_expression_trees':len(ff)})
H=(D/'calendar.hpp').read_text();native_months=re.findall(r"'([A-Z]{3})'",S[:S.index('contains')]);cpp_months=re.findall(r'"([A-Z]{3})"',H);assert native_months==cpp_months
assert 'std::array<int, 12> days_per_month{31, 28, 31, 30, 31, 30,' in H and '31, 31, 30, 31, 30, 31};' in H
assert 'double seconds_since_previous_date' in C and 'double& seconds_since_current_date' in C
assert not re.search(r'\bthread_local\b|f2k::|KOKKOS_LAMBDA|\bfma\s*\(',C)
result={'scope':'All five Native calendar routines and two immutable tables, full Native/JAX read before outline.','canonical_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [D/'calendar.cpp',D/'calendar.hpp']},'reviews':reviews,'source_notes':['Exact Native truncating INTEGER expressions and source assignment order.','Elapsed FLOOR precedes integer conversion; seconds use Native time_precision double, independent of core kind.','Native repeated leap predicates, one-based prefix bounds, upper-day fatal branch/message retained; JAX omitted fatal not imported.','Only immutable Native tables and call-local scalars. No new state owner or lower bridge.'],'provenance':'Owner source review; independent review and separate cleanup tracked separately.','runtime_status':'Draft only. No builds/tests/CPU/GPU/BFB claim; all inputs, rollovers, invalid-day error branch await integration.'};(D/'SOURCE_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print('Calendar source audit passed:',len(reviews),'routines.')
