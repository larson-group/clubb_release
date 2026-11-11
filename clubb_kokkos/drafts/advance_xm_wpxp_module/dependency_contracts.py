"""Record actual public Native lower contracts without editing the f2k owner."""
from pathlib import Path
from dataclasses import asdict
import json,re,sys
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import find_procedure
names=['calc_wp3_on_wp2','calc_turb_adv_range','diffusion_zm_lhs','term_ma_zm_lhs',
       'diffusion_zt_lhs','term_ma_zt_lhs','xpyp_term_ta_pdf_rhs','xpyp_term_ta_pdf_lhs',
       'xpyp_term_ta_pdf_lhs_godunov','monotonic_turbulent_flux_limit','pos_definite_adj',
       'fill_holes_vertical_api','clip_covar','stats_update_budget_2d']
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
 arrays=[a for a in p.arguments if a.dims];types=[]
 for a in p.arguments:
  if a.dims:typ='const Array'+str(arrays.index(a))+'&'
  elif a.kind in {'grid','implicit_coefs_terms','nu_vertical_res_dep','err_info_type'}:typ=('const ' if a.intent=='in' else '')+'clubb_kokkos::'+a.kind+'&'
  elif a.kind=='stats_type':typ='stats_type&'
  elif a.kind=='character':typ='const std::string&'
  else:typ={'integer':'int','logical':'bool','real':'clubb_kokkos::real_t','double':'double'}[a.kind]+('&' if a.intent!='in' else '')
  assert not a.optional
  types.append(typ+' '+a.name)
 head+='// '+p.source+':'+p.name+'; all arguments required.\n'
 head+='template<'+','.join('class Array'+str(i) for i in range(len(arrays)))+'>\n'
 head+='void '+record['cpp_name']+'('+',\n'.join(types)+');\n'
head+='} // namespace f2k\n'
(folder/'f2k_dependencies.hpp').write_text(head)
(folder/'lower_contracts.json').write_text(json.dumps(records,indent=2)+'\n')
print('Missing public bridge definitions:',', '.join(r['cpp_name'] for r in records if not r['existing_bridge']))
