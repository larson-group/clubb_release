"""Record actual public Native lower contracts without editing the f2k owner."""
from pathlib import Path
from dataclasses import asdict
import json,re,sys
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import find_procedure
names=['calc_wp3_on_wp2','diffusion_zm_lhs','term_ma_zm_lhs','xpyp_term_ta_pdf_rhs','xpyp_term_ta_pdf_lhs','xpyp_term_ta_pdf_lhs_godunov','xpyp_term_ta_pdf_rhs_godunov','fill_holes_vertical_api','clip_covar','clip_variance','stats_update_budget_2d']
existing=(root/'clubb_kokkos/src/f2k/procedure_bridge.hpp').read_text()
records=[];head='''// Declaration-only Native lower boundaries for this isolated draft.
// Root must add missing definitions through the existing f2k generator owner,
// then remove this staging header. These declarations are not bridge implementations.
#pragma once
#include "f2k/procedure_bridge.hpp"
namespace f2k {
'''
for name in names:
 p=find_procedure(name);record=asdict(p);record['cpp_name']='stats_update_budget' if name=='stats_update_budget_2d' else name
 record['existing_bridge']=bool(re.search(r'inline void '+record['cpp_name']+r'\s*\(',existing))
 records.append(record)
 if record['existing_bridge']:continue
 arrays=[a for a in p.arguments if a.dims and not a.optional];types=[]
 for a in p.arguments:
  if a.dims:typ='const clubb_kokkos::Real'+str(len(a.dims))+'D*' if a.optional else 'const Array'+str(arrays.index(a))+'&'
  elif a.kind in {'grid','implicit_coefs_terms','nu_vertical_res_dep','err_info_type'}:typ=('const ' if a.intent=='in' else '')+'clubb_kokkos::'+a.kind+'&'
  elif a.kind=='stats_type':typ='stats_type&'
  elif a.kind=='character':typ='const std::string&'
  else:typ={'integer':'int','logical':'bool','real':'clubb_kokkos::real_t','double':'double'}[a.kind]+('&' if a.intent!='in' else '')
  if a.optional and not a.dims:typ='const '+typ+'*'
  types.append(typ+' '+a.name+(' = nullptr' if a.optional else ''))
 head+='// '+p.source+':'+p.name+('; Native optional threshold_hi retains presence semantics.\n' if any(a.optional for a in p.arguments) else '; all arguments required.\n')
 head+='template<'+','.join('class Array'+str(i) for i in range(len(arrays)))+'>\n'
 head+='void '+record['cpp_name']+'('+',\n'.join(types)+');\n'
extra=find_procedure('stats_finalize_budget_2d');record=asdict(extra);record['cpp_name']='stats_finalize_budget';record['existing_bridge']='Only the3-argument default exists; explicit fourth optional bool remains pending';records.append(record)
head+='// Native stats_finalize_budget_2d optional l_count_sample; explicit-value overload.\n'
head+='template<class Array0> void stats_finalize_budget(const std::string& name,const Array0& values,stats_type& stats,bool l_count_sample);\n'
head+='} // namespace f2k\n'
(folder/'f2k_dependencies.hpp').write_text(head)
(folder/'lower_contracts.json').write_text(json.dumps(records,indent=2)+'\n')
print('Missing public bridge definitions:',', '.join(r['cpp_name'] for r in records if not r['existing_bridge']))
