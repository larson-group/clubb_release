from pathlib import Path
import sys,re
root=Path('/home/guntherhuebler/Code/clubb');stage=Path('/tmp/clubb_adg1_next')
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import find_procedure
source=(root/'src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90').read_text()
names=['ADG1_pdf_driver','ADG1_ADG2_responder_params']
blocks={n:re.search(r'^  subroutine '+n+r'\(.*?^  end subroutine '+n+r'\b',source,re.M|re.S).group() for n in names}
procs={n:find_procedure(n) for n in names}
def ctype(a):
    if a.kind=='real':return 'Real'+str(len(a.dims))+'D' if a.dims else 'real_t'
    if a.kind=='logical':return 'bool'
    if a.kind=='integer':return 'int'
    if a.kind=='err_info_type':return 'err_info_type&'
    raise ValueError(a)
def signature(n):
    args={a.name.lower():a for a in procs[n].arguments}
    head='\n'.join(l.split('!')[0] for l in blocks[n].splitlines());head=head[:head.index(')')]
    groups=[]
    for line in head.splitlines():
        code=re.sub(r'^\s*subroutine\s+\w+\(','',line.replace('&',''))
        fields=[x.strip().lower() for x in code.split(',') if x.strip()]
        if fields:groups.append((', '.join(ctype(args[x])+' '+args[x].name for x in fields),'+'.join(dict.fromkeys(args[x].intent for x in fields))))
    labels={'in':'In','out':'Out','inout':'InOut'}
    return 'void '+n+'('+''.join(g+(', // '+ '/'.join(labels[j] for j in intents.split('+'))+'\n    ' if i<len(groups)-1 else ')') for i,(g,intents) in enumerate(groups))
def call(callee,index=0,prefix=''):
    b=blocks['ADG1_pdf_driver'];lines=b.splitlines();starts=[i for i,l in enumerate(lines) if re.search(r'\bcall\s+'+callee+r'\(',l,re.I)]
    i=starts[index];out=[];depth=0;first=True
    while i<len(lines):
        line=lines[i];code=line.split('!')[0];depth+=code.count('(')-code.count(')')
        code=re.sub(r'^\s*call\s+','',code,flags=re.I).strip().replace('&','').strip()
        code=re.sub(r'(\w+)\s*\(\s*:\s*,\s*:\s*,\s*sclr\s*\)',r'Kokkos::subview(\1, Kokkos::ALL, Kokkos::ALL, sclr)',code)
        code=code.replace('sclr_tol(sclr)','sclr_tol_host(sclr)')
        if first:code=prefix+code;first=False
        comment=line.split('!',1)[1].strip() if '!' in line else ''
        comment=comment.replace('In/Out','InOut')
        if depth==0:code+=';'
        out.append('  '+code+(' // '+comment if comment else ''))
        i+=1
        if depth==0:break
    return '\n'.join(out)
header='''// Compute ADG1 PDF component means and variances from grid-box moments.
//
// ADG1_pdf_driver first closes vertical velocity, then computes responder
// parameters for total water, liquid-water potential temperature, horizontal
// winds and optional passive scalars. Component means retain each field's
// units; variances have squared units and alpha/mixture factors are dimensionless.
// All profile Views use execution-space, column-first LayoutLeft storage;
// scalar fields add a final species dimension. The caller owns outputs and
// err_info. Vertical-velocity closure remains an explicit temporary f2k call.
// Source: src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90:ADG1_pdf_driver.
#pragma once

#include "array_types.hpp"
#include "err_info_type_module.hpp"

namespace clubb_kokkos {

'''+signature('ADG1_pdf_driver')+';\n\n} // namespace clubb_kokkos\n'
(stage/'adg1_adg2_3d_luhar_pdf.hpp').write_text(header)
cpp='''// ADG1 closure orchestration and responder parameters; see the public header.
#include "adg1_adg2_3d_luhar_pdf.hpp"
#include "constants_clubb.hpp"
#include "error_code.hpp"
#include "f2k/procedure_bridge.hpp"
#include <iostream>

namespace clubb_kokkos {

namespace {

'''+signature('ADG1_ADG2_responder_params')+';\n\n} // namespace\n\n'
body='''  // Calculate ADG1 component parameters from column/level moments, using
  // the vertical-velocity closure followed by rt, thl, u, v and scalar
  // responders. Return component means/variances and update native error state.

  // Local Variables
  Real2D w_1_n("ADG1_pdf_driver::w_1_n", ngrdcol, nz); // Normalized component1 w mean [-]
  Real2D w_2_n("ADG1_pdf_driver::w_2_n", ngrdcol, nz); // Normalized component2 w mean [-]

  // Calculate mixture fraction and PDF component means and variances of w.
  // TODO: Replace this native lower call when ADG1_w_closure is ported.
'''+call('ADG1_w_closure',prefix='f2k::')+'''

  if (f2k::clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      return;
    }
  }

  // Calculate the PDF component means and variances of rt.
'''+call('ADG1_ADG2_responder_params',0)+'''

  // Calculate the PDF component means and variances of thl.
'''+call('ADG1_ADG2_responder_params',1)+'''

  // Calculate the PDF component means and variances of u wind.
'''+call('ADG1_ADG2_responder_params',2)+'''

  // Calculate the PDF component means and variances of v wind.
'''+call('ADG1_ADG2_responder_params',3)+'''

  // Calculate the PDF component means and variances of passive scalars.
  if (l_scalar_calc) {
    // Per-species tolerances cross from the execution space only at this
    // host dispatch boundary; scalar fields remain column-first device Views.
    const auto sclr_tol_host = to_host_space(sclr_tol);
    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
'''+call('ADG1_ADG2_responder_params',4)+'''
    }
  }
'''
cpp+=signature('ADG1_pdf_driver')+' {\n'+body+'}\n\nnamespace {\n\n'+signature('ADG1_ADG2_responder_params')+' {\n'
# Preserve the complete native scientific description and its historical notes.
b=blocks['ADG1_ADG2_responder_params'];start=b.index('    ! Description:');end=b.index('    use constants_clubb')
cpp+='\n'.join('  //'+l.strip()[1:] if l.strip().startswith('!') else '' for l in b[start:end].splitlines())+'\n'
cpp+='''  //-------------------------- Begin Code --------------------------
  if (f2k::clubb_at_least_debug_level_api(2)) {
    // Integer ANY reduction leaves valid data on the execution space.
    int below_tolerance = 0;
    Kokkos::parallel_reduce("ADG1_ADG2_responder_params::variance_check", ngrdcol * nz,
      KOKKOS_LAMBDA(const int cell, int& count) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        if (xp2(i,k) < x_tol * x_tol) ++count;
      }, below_tolerance);
    if (below_tolerance) {
      std::cerr << err_info.err_header_global << '\\n';
      std::cerr << "Error in ADG1_ADG2_responder_params: xp2 must not be less "
                   "than its tolerance value.\\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }

  Kokkos::parallel_for("ADG1_ADG2_responder_params::component_parameters", ngrdcol * nz,
    KOKKOS_LAMBDA(const int cell) {
      const int i = cell % ngrdcol;
      const int k = cell / ngrdcol;
      real_t width_factor_1; // Width factor relating to PDF component1 [-]

      // This if-statement has been commented out in favor of the error check above.
      // These fields should never be below tolerance, and these tolerances are too large
      // to make this if-case a proper asymptote of the computation in the else-case.
      // See https://github.com/larson-group/clubb/pull/1276 for more detail.
      // if (wp2(i,k) < w_tol_sqd || xp2(i,k) < x_tol*x_tol) {
      //   // Vertical velocity doesn't vary. Variable x is treated as a single
      //   // Gaussian in this situation.
      //   x_1(i,k) = xm(i,k);
      //   x_2(i,k) = xm(i,k);
      //   varnce_x_1(i,k) = xp2(i,k);
      //   varnce_x_2(i,k) = xp2(i,k);
      //   alpha_x(i,k) = one_half;
      // } else {

      x_1(i,k) = xm(i,k) - wpxp(i,k) / (sqrt_wp2(i,k) * w_2_n(i,k));
      x_2(i,k) = xm(i,k) - wpxp(i,k) / (sqrt_wp2(i,k) * w_1_n(i,k));

      alpha_x(i,k) = one_half * (one - wpxp(i,k) * wpxp(i,k)
          / ((one - sigma_sqd_w(i,k)) * wp2(i,k) * xp2(i,k)));

      // Native MIN/MAX select their second argument on an equality tie.
      alpha_x(i,k) = alpha_x(i,k) < one ? alpha_x(i,k) : one;
      alpha_x(i,k) = alpha_x(i,k) > zero_threshold ? alpha_x(i,k) : zero_threshold;

      width_factor_1 = two_thirds * beta(i)
          + two * mixt_frac(i,k) * (one - two_thirds * beta(i));

      // Vince Larson multiplied original expressions by width_factor_1,2
      // to generalize scalar skewnesses. 05Nov2003.
      varnce_x_1(i,k) = width_factor_1 * xp2(i,k) * alpha_x(i,k) / mixt_frac(i,k);
      varnce_x_2(i,k) = (two - width_factor_1) * xp2(i,k)
          * alpha_x(i,k) / (one - mixt_frac(i,k));

      // } // Former single-Gaussian branch remains inactive.
    });
}

} // namespace

} // namespace clubb_kokkos
'''
(stage/'adg1_adg2_3d_luhar_pdf.cpp').write_text(cpp)
