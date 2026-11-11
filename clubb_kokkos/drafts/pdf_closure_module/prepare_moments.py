from pathlib import Path
import re,json,sys
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import find_procedure,statements
from generate_derived_types import split_top
source=(root/'src/CLUBB_core/pdf_closure_module.F90').read_text()
names=['calc_wp4_pdf','calc_wp2xp_pdf','calc_wpxp2_pdf','calc_wpxpyp_pdf']
blocks={n:re.search(r'^  subroutine '+n+r'\(.*?^  end subroutine '+n+r'\b',source,re.M|re.S).group() for n in names}
outline={'authority':'src/CLUBB_core/pdf_closure_module.F90','routines_in_native_order':names,'contracts':{n:{'arguments':[(a.name,a.kind,a.intent,a.dims) for a in find_procedure(n).arguments],'options':[],'branches':[],'lower_calls':[],'work_mapping':'independent column/level cells','outputs':'last required array overwritten in every cell'} for n in names},'comments':'Retain full Native equation/parameter descriptions, argument units and calculation comment.','adaptations':['Column-first LayoutLeft Kokkos loops replace Native nested loops.','Integer squares expand to products, fourth powers to squared squares; no new guards/temporaries/FMA.','At promotion public definitions must interrupt/reopen anonymous private-helper scopes to retain Native definition order.'],'runtime':'not compiled/integrated/tested'}
(folder/'SOURCE_OUTLINE.json').write_text(json.dumps(outline,indent=2)+'\n')
def powers(s):
 while '**' in s:
  op=s.index('**');end=op;j=op-1
  while s[j].isspace():j-=1
  if s[j]==')':
   depth=1;j-=1
   while depth:
    depth+=(s[j]==')')-(s[j]=='(');j-=1
   start=j+1
   while j>=0 and (s[j].isalnum() or s[j]=='_'):j-=1;start=j+1
  else:
   while j>=0 and (s[j].isalnum() or s[j]=='_'):j-=1
   start=j+1
  m=re.match(r'\*\*\s*([24])\b',s[op:]);assert m,s[op:]
  base=s[start:end].strip();square=f'({base} * {base})'
  replacement=square if m[1]=='2' else f'({square} * {square})'
  s=s[:start]+replacement+s[op+m.end():]
 return s
cpp='''// Higher moments integrated over the two-component normal PDF.
#include "pdf_closure_module.hpp"
#include "constants_clubb.hpp"

namespace clubb_kokkos {

''';hpp='''// Integrate fourth/third moments over the two-component normal PDF.
// Source-only draft: calc_wp4_pdf, calc_wp2xp_pdf, calc_wpxp2_pdf and calc_wpxpyp_pdf.
// Column-first View arguments retain Native field units; not compiled or tested.
#pragma once
#include "array_types.hpp"

namespace clubb_kokkos {

'''
for n in names:
 native=blocks[n];proc=find_procedure(n);out=proc.arguments[-1].name
 sig=native[:native.index(')')+1];rows=[]
 for l in sig.splitlines():
  l=l.split('!')[0].strip().replace('&','')
  l=re.sub(r'^subroutine\s+'+n+r'\(\s*','',l).rstrip(' )')
  args=[x.strip() for x in l.split(',') if x.strip()]
  rows.append(', '.join(('int ' if a.kind=='integer' else 'Real2D ')+name for name in args for a in proc.arguments if a.name==name))
 signature='void '+n+'('+('\n'.join(row+(',' if i<len(rows)-1 else ')')+(' // In' if i<len(rows)-1 else '') for i,row in enumerate(rows)))
 hpp+=signature+';\n\n';cpp+=signature+' {\n'
 comments=native[native.index('    ! Description:'):native.index('    use constants_clubb')]
 cpp+='\n'.join('  //'+l.split('!',1)[1] if '!' in l else '' for l in comments.splitlines())+'\n'
 cpp+='  // Arrays are indexed (column, level).\n'
 for l in native[native.index('    ! Input Variables'):native.index('    ! Local Variables')].splitlines():
  m=re.match(r'\s*(\w+)\s*,?\s*&?\s*!\s*(.*)',l)
  if m:cpp+='  // '+m[1]+': '+re.sub(r' {2,}', ' ', m[2])+'\n'
 assignment=next(s for s in statements(native) if re.match(r'\s*'+out+r'\s*\(',s) and '=' in s)
 assignment=powers(assignment);assignment=re.sub(r'\bsqrt\s*\(','Kokkos::sqrt(',assignment)
 assignment=re.sub(r'\(\s*i\s*,\s*k\s*\)','(i, k)',assignment)
 calc=next(l.split('!',1)[1].strip() for l in native.splitlines() if '! Calculate ' in l)
 cpp+='\n  // '+calc+'\n'
 cpp+='  Kokkos::parallel_for(\n      "'+n+'", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {\n        const int i = cell % ngrdcol;\n        const int k = cell / ngrdcol;\n        '+assignment+';\n      });\n}\n\n'
cpp+='} // namespace clubb_kokkos\n';hpp+='} // namespace clubb_kokkos\n'
(folder/'pdf_closure_module.cpp').write_text(cpp);(folder/'pdf_closure_module.hpp').write_text(hpp)
print('Drafted',names)
