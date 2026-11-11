from pathlib import Path
import sys,re,json,ast
root=Path('/home/guntherhuebler/Code/clubb');stage=Path('/tmp/clubb_adg1_next')
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import find_procedure,statements
from generate_derived_types import split_top
source=(root/'src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90').read_text()
names=['ADG1_pdf_driver','ADG1_ADG2_responder_params']
blocks={n:re.search(r'^  subroutine '+n+r'\(.*?^  end subroutine '+n+r'\b',source,re.M|re.S).group() for n in names}
cpp=(stage/'adg1_adg2_3d_luhar_pdf.cpp').read_text();header=(stage/'adg1_adg2_3d_luhar_pdf.hpp').read_text()
clean=lambda s:re.sub(r'//[^\n]*|/\*.*?\*/','',s,flags=re.S)
code=clean(cpp)
defs=re.findall(r'void\s+(\w+)\s*\(([^)]*)\)\s*\{',code)
assert [n for n,a in defs]==names
for n,args in defs:
    assert [x.strip().split()[-1] for x in split_top(args)]==[a.name for a in find_procedure(n).arguments]
ha=re.search(r'void\s+ADG1_pdf_driver\s*\(([^)]*)\)',clean(header)).group(1)
assert [x.strip().split()[-1] for x in split_top(ha)]==[a.name for a in find_procedure(names[0]).arguments]
def calls(s,cpp=False):
    s=clean(s) if cpp else '\n'.join(l.split('!')[0] for l in s.splitlines())
    pat=r'(?<![\w:])(?:f2k::)?(ADG1_w_closure|ADG1_ADG2_responder_params)\s*\('
    if not cpp:pat=r'call\s+'+pat
    result=[]
    for m in re.finditer(pat,s):
        start=m.end();pos=start;depth=1
        while depth:
            depth+=(s[pos]=='(')-(s[pos]==')');pos+=1
        args=s[start:pos-1]
        if cpp and re.search(r'\b(?:int|Real\dD|const)\b',args):continue
        result.append((m[1],split_top(args.replace('&',''))))
    return result
nc=calls(blocks[names[0]]);cc=calls(cpp,True)
assert [n for n,a in nc]==[n for n,a in cc]
def norm(a,cpp=False):
    if cpp:
        a=re.sub(r'Kokkos::subview\(\s*(\w+)\s*,\s*Kokkos::ALL\s*,\s*Kokkos::ALL\s*,\s*sclr\s*\)',r'\1(:,:,sclr)',a)
        a=a.replace('sclr_tol_host','sclr_tol')
    return re.sub(r'\s','',a).lower()
for (nn,na),(cn,ca) in zip(nc,cc):assert [norm(a) for a in na]==[norm(a,True) for a in ca],nn
# Compare the complete formula trees of every unchanged responder assignment.
# Clamp is the separately reviewed native MIN/MAX language adaptation.
fn=blocks[names[1]]
assignments={}
for st in statements(fn):
    if re.match(r'(?:x_[12]|alpha_x|width_factor_1|varnce_x_[12])(?:\(i,k\))?\s*=',st):
        lhs,rhs=st.split('=',1)
        if 'max(' not in rhs:assignments[lhs.strip()]=rhs
for lhs,rhs in assignments.items():
    lhs_rx=re.escape(lhs).replace(',',r',\s*')
    target=re.search(lhs_rx+r'\s*=\s*(.*?);',code,re.S).group(1)
    a=ast.dump(ast.parse(rhs.strip(),mode='eval'),include_attributes=False)
    b=ast.dump(ast.parse(' '.join(target.split()),mode='eval'),include_attributes=False)
    assert a==b,(lhs,rhs,target)
assert 'rt_tol' in cc[1][1] and 'thl_tol' in cc[2][1] and 'thl_tol' in cc[3][1] and 'thl_tol' in cc[4][1]
assert cpp.count('if (any_equal(err_info.err_code, clubb_fatal_error))')==1
assert 'Kokkos::deep_copy(err_info.err_code, clubb_fatal_error)' in code
assert cpp.index('if (xp2(i, k) < x_tol * x_tol)') < cpp.index('"ADG1_ADG2_responder_params::component_parameters"')
for p in ['adg1_adg2_3d_luhar_pdf.cpp','adg1_adg2_3d_luhar_pdf.hpp']:
    text=(stage/p).read_text();s=clean(text);s=re.sub(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'','',s)
    for a,b in [('(',')'),('{','}'),('[',']')]:
        depth=0
        for char in s:
            depth+=(char==a)-(char==b);assert depth>=0,(p,a)
        assert depth==0,(p,a,depth)
    assert all(len(l)<=100 and l.rstrip()==l and '\t' not in l for l in text.splitlines())
result={'source':'src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90:32–240,1072–1233','routines':names,'arguments':{n:len(find_procedure(n).arguments) for n in names},'source_definition_order':True,'calls_in_source_order':len(nc),'every_call_argument_mapping':True,'responder_expression_trees_compared':len(assignments),'clamp_manual_review':'min(alpha,one), then max(result,zero_threshold), retaining second operand on equal tie','manual_review':['source signatures, shapes/intents/required arguments, call grouping/order','sole parent fatal gate immediately after w closure, no new responder gates','all-column debug2 variance check and untouched outputs on return','rt/thl/u/v field mappings and native thl_tol for both winds','explicit l_scalar_calc branch, per-species tolerance host boundary, untouched outputs when skipped','source equations/denominator order, scalar width temporary per GPU work item','scientific derivation,2003 width-factor notes and inactive branch/PR1276 rationale','only w closure remains lower f2k boundary; no stats/global state/new UDT','comment/spacing/formatting cleanup and natural-language adaptations'], 'builds_or_numerical_tests':'not run; root owns later integration/model comparison'}
(stage/'source_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
