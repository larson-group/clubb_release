"""Prepare the isolated PDF utility draft from current Native contracts/comments.

Only files in this script's draft directory are written. This prepares source;
it does not compile, execute or establish numerical equivalence.
"""
from pathlib import Path
import re
import json
import hashlib

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
SOURCE = ROOT / "src/CLUBB_core/pdf_utilities.F90"
NATIVE = SOURCE.read_text()
ROUTINES = {}
for match in re.finditer(r"^  (elemental function|subroutine) (\w+)\b", NATIVE, re.M):
    name = match.group(2)
    kind = "function" if match.group(1).startswith("elemental") else "subroutine"
    end = NATIVE.index("end " + kind + " " + name, match.start())
    end = NATIVE.index("\n", end)
    block = NATIVE[match.start():end]
    signature = ""
    for line in block.splitlines():
        code = line.split("!", 1)[0].strip()
        signature += " " + code.strip("& ")
        if not code.endswith("&"):
            break
    args = [a.strip() for a in re.search(r"\b" + name + r"\s*\((.*?)\)", signature).group(1).split(",")]
    ROUTINES[name] = {"block": block, "args": args, "line": NATIVE.count("\n", 0, match.start()) + 1, "kind": kind}

RESULTS = dict(zip("mean_L2N mean_L2N_dp stdev_L2N stdev_L2N_dp corr_NL2NN corr_NL2NN_dp corr_LL2NN corr_LL2NN_dp compute_mean_binormal compute_variance_binormal calc_corr_chi_x calc_corr_eta_x calc_corr_rt_x calc_corr_thl_x calc_xp2".split(), "mu_x_n mu_x_n sigma_x_n sigma_x_n corr_x_y_n corr_x_y_n corr_x_y_n corr_x_y_n xm xp2 corr_chi_x_i corr_eta_x_i corr_rt_x_i corr_thl_x_i xp2".split()))
BODY_START = {}
for name in RESULTS:
    if name.startswith("mean_L2N"):
        BODY_START[name] = "    ! Find the mean of ln x"
    elif name.startswith("stdev_L2N"):
        BODY_START[name] = "    ! Find the standard deviation"
    elif name.startswith("corr_"):
        BODY_START[name] = "    !-------------------------------------- Begin Code"
    elif name.startswith("compute_"):
        BODY_START[name] = "    !----- Begin Code"
    elif name.startswith("calc_corr"):
        BODY_START[name] = "    ! Calculate the correlation"
    else:
        BODY_START[name] = "    ! Calculate overall variance"


def description(block):
    """Retain the Native procedure description, equations and references."""
    if "! Description:" not in block:
        return []
    part = block[block.index("! Description:"):]
    part = part[:part.index("\n    use ")]
    return ["  //" + line.split("!", 1)[1] if "!" in line else "" for line in part.splitlines()]


def arguments(info):
    """Retain Native argument meanings/units beside the C++ implementation."""
    comments = {}
    for line in info["block"].splitlines():
        if "!" not in line:
            continue
        code, comment = line.split("!", 1)
        match = re.match(r"\s*(\w+)\s*[, &]*$", code)
        if match and match.group(1) in info["args"]:
            comments[match.group(1)] = comment.strip()
    return ["  // " + name + ": " + comments[name] for name in info["args"] if name in comments]


def expression(code):
    """Translate only the scalar syntax present in these Native procedures."""
    code = re.sub(r"real\(\s*max_mag_correlation,\s*kind\s*=\s*dp\s*\)", "static_cast<double>(max_mag_correlation)", code, flags=re.I)
    code = re.sub(r"\(\s*(mu_x_[12])\s*-\s*xm\s*\)\s*\*\*\s*2", lambda m: "((" + m.group(1) + " - xm) * (" + m.group(1) + " - xm))", code)
    code = re.sub(r"\b([a-zA-Z_]\w*)\s*\*\*\s*2\b", lambda m: "(" + m.group(1) + " * " + m.group(1) + ")", code)
    if "**" in code:
        raise ValueError("Untranslated power: " + code)
    code = re.sub(r"tiny\(\s*mu_x\s*\)", "std::numeric_limits<real_t>::min()", code, flags=re.I)
    code = re.sub(r"\b(log|sqrt|exp|abs|max)\s*\(", lambda m: "Kokkos::" + m.group(1) + "(", code)
    return code.replace(".and.", "&&").replace(".or.", "||")


def scalar_body(name, info):
    """Keep scalar assignment and branch order while replacing Fortran syntax."""
    block = info["block"]
    part = block[block.index(BODY_START[name]):]
    part = part[:part.index("\n  end function")]
    lines = []
    pending = ""
    indent = 2
    for raw in part.splitlines():
        if raw.lstrip().startswith("!"):
            if not raw.lstrip().startswith("!$"):
                lines.append(" " * indent + "//" + raw.split("!", 1)[1])
            continue
        code = raw.split("!", 1)[0].strip()
        if not code:
            if not pending and (not lines or lines[-1] != ""):
                lines.append("")
            continue
        pending += " " + code.strip("& ")
        if code.endswith("&"):
            continue
        code = pending.strip()
        pending = ""
        if code.lower() == "return":
            lines.append(" " * indent + "return " + RESULTS[name] + ";")
        elif re.match(r"endif|end if", code, re.I):
            indent -= 2
            lines.append(" " * indent + "}")
        elif re.match(r"elseif", code, re.I):
            indent -= 2
            cond = re.sub(r"^elseif\s*\((.*)\)\s*then$", r"\1", code, flags=re.I)
            lines.append(" " * indent + "} else if (" + expression(cond) + ") {")
            indent += 2
        elif code.lower() == "else":
            indent -= 2
            lines.append(" " * indent + "} else {")
            indent += 2
        elif re.match(r"if\s*\(", code, re.I):
            cond = re.sub(r"^if\s*\((.*)\)\s*then$", r"\1", code, flags=re.I)
            lines.append(" " * indent + "if (" + expression(cond) + ") {")
            indent += 2
        else:
            lines.append(" " * indent + expression(code) + ";")
    if pending or indent != 2:
        raise ValueError((name, pending, indent))
    return lines


HEADER = '''// Convert PDF component moments and correlations between normal/lognormal variables.
//
// Lognormal conversions: mean_L2N/mean_L2N_dp and stdev_L2N/stdev_L2N_dp give
// moments of ln(x); corr_NL2NN/corr_NL2NN_dp and corr_LL2NN/corr_LL2NN_dp map
// correlations into normal space; corr_NN2NL and corr_NN2LL perform the inverse.
// Binormal moments: compute_mean_binormal, compute_variance_binormal and
// calc_comp_corrs_binormal combine the two normal PDF components.
// Coordinate conversions: calc_corr_chi_x/calc_corr_eta_x and their inverses
// calc_corr_rt_x/calc_corr_thl_x map correlations between thermodynamic and PDF
// coordinates. calc_xp2 combines lognormal/zero component variances;
// smooth_corr_quotient regularizes a correlation denominator smoothly.
//
// This mirrors all 19 public routines in src/CLUBB_core/pdf_utilities.F90.
// Scalar elemental forms are device callable; native array procedures use
// execution-space Views indexed (column, level), with LayoutLeft. The three
// rank-two elemental overloads at the end serve existing PDF closure callers.
// Inputs are required and conformable, as in Native. Storage belongs to each
// call. Distinct _dp behavior is retained even when real_t is already double.
// This isolated draft has not been compiled, integrated or numerically tested.
// Scientific references are retained beside the corresponding procedures.
#pragma once

#include "array_types.hpp"
#include "constants_clubb.hpp"
#include "advance_helper_module.hpp"
#include <limits>

namespace clubb_kokkos {

// Forward declaration permits the native calc_comp_corrs_binormal call order.
inline void smooth_corr_quotient(int ngrdcol, int nz, Real2D numerator, Real2D denominator,
                                 real_t denom_thresh, Real2D quotient);
'''

ARRAYS = {}
ARRAYS['corr_NN2NL'] = r'''
inline void corr_NN2NL(int nz, int ngrdcol,
                       Real2D corr_x_y_n, Real2D sigma_y_n, Real2D y_sigma2_on_mu2, // In
                       Real2D corr_x_y) {                                        // Out
{DESCRIPTION}

  // Input Variables
{ARGUMENTS}

  //-------------------------------------- Begin Code --------------------------------------
  // Find the correlation of x and y for the ith component of the PDF.
  // When sigma_y=0 and mu_y>0, y_sigma2_on_mu2=0 and sigma_y_n=0.
  // The correlations are then undefined; preserve the native limiting value.
  Kokkos::parallel_for("corr_NN2NL", ngrdcol*nz, KOKKOS_LAMBDA(const int cell) {
    const int j = cell % ngrdcol;
    const int k = cell / ngrdcol;
    if (sigma_y_n(j,k) > zero) {
      // This ratio is already positive here. The tiny floor prevents an
      // underflowed denominator from causing division by zero.
      corr_x_y(j,k) = corr_x_y_n(j,k) * sigma_y_n(j,k)
                       / Kokkos::sqrt(Kokkos::max(y_sigma2_on_mu2(j,k),
                                                 std::numeric_limits<real_t>::min()));
    } else {
      // sqrt(log(1+r)/r) tends to1 as r approaches0 from the right.
      corr_x_y(j,k) = corr_x_y_n(j,k);
    }
    // Preserve the native ordered upper/lower clipping checks.
    if (corr_x_y(j,k) > max_mag_correlation) {
      corr_x_y(j,k) = max_mag_correlation;
    } else if (corr_x_y(j,k) < -max_mag_correlation) {
      corr_x_y(j,k) = -max_mag_correlation;
    }
  });
}
'''
ARRAYS['corr_NN2LL'] = r'''
inline void corr_NN2LL(int nz, int ngrdcol,
                       Real2D corr_x_y_n, Real2D sigma_x_n, Real2D sigma_y_n, // In
                       Real2D x_sigma2_on_mu2, Real2D y_sigma2_on_mu2,       // In
                       Real2D corr_x_y) {                                  // Out
{DESCRIPTION}

  // Input Variables
{ARGUMENTS}

  //-------------------------------------- Begin Code --------------------------------------
  // Either zero transformed standard deviation makes both correlations
  // undefined. Return the supplied correlation in that native limiting case.
  Kokkos::parallel_for("corr_NN2LL", ngrdcol*nz, KOKKOS_LAMBDA(const int cell) {
    const int j = cell % ngrdcol;
    const int k = cell / ngrdcol;
    if (sigma_x_n(j,k) > zero && sigma_y_n(j,k) > zero) {
      // Native alternate form:
      // (exp(sigma_x_n*sigma_y_n*corr_x_y_n)-one)
      //   /(sqrt(exp(sigma_x_n**2)-one)*sqrt(exp(sigma_y_n**2)-one)).
      corr_x_y(j,k) = (Kokkos::exp(sigma_x_n(j,k)*sigma_y_n(j,k)*corr_x_y_n(j,k))-one)
                       / Kokkos::sqrt(x_sigma2_on_mu2(j,k)*y_sigma2_on_mu2(j,k));
    } else {
      corr_x_y(j,k) = corr_x_y_n(j,k);
    }
    // Clip inconsistent component moments to a realizable correlation.
    if (corr_x_y(j,k) > max_mag_correlation) {
      corr_x_y(j,k) = max_mag_correlation;
    } else if (corr_x_y(j,k) < -max_mag_correlation) {
      corr_x_y(j,k) = -max_mag_correlation;
    }
  });
}
'''
ARRAYS['calc_comp_corrs_binormal'] = r'''
inline void calc_comp_corrs_binormal(int nz, int ngrdcol, // In
                                     Real2D xpyp, Real2D xm, Real2D ym, // In
                                     Real2D mu_x_1, Real2D mu_x_2,      // In
                                     Real2D mu_y_1, Real2D mu_y_2,      // In
                                     Real2D sigma_x_1_sqd,              // In
                                     Real2D sigma_x_2_sqd,              // In
                                     Real2D sigma_y_1_sqd,              // In
                                     Real2D sigma_y_2_sqd,              // In
                                     Real2D mixt_frac,                  // In
                                     Real2D corr_x_y_1,                 // Out
                                     Real2D corr_x_y_2) {               // Out
{DESCRIPTION}

  // Input/Output Variables
{ARGUMENTS}

  //------------------------ Local Variables ------------------------
  Real2D numerator("calc_comp_corrs_binormal::numerator", ngrdcol,nz);
  Real2D denominator("calc_comp_corrs_binormal::denominator", ngrdcol,nz);

  //-------------------------------------- Begin Code --------------------------------------
  Kokkos::parallel_for("calc_comp_corrs_binormal::quotient_terms", ngrdcol*nz,
    KOKKOS_LAMBDA(const int cell) {
      const int i = cell % ngrdcol;
      const int k = cell / ngrdcol;
      numerator(i,k) = xpyp(i,k)
        - mixt_frac(i,k)*(mu_x_1(i,k)-xm(i,k))*(mu_y_1(i,k)-ym(i,k))
        - (one-mixt_frac(i,k))*(mu_x_2(i,k)-xm(i,k))*(mu_y_2(i,k)-ym(i,k));
      denominator(i,k) = mixt_frac(i,k)*Kokkos::sqrt(sigma_x_1_sqd(i,k)*sigma_y_1_sqd(i,k))
        + (one-mixt_frac(i,k))*Kokkos::sqrt(sigma_x_2_sqd(i,k)*sigma_y_2_sqd(i,k));
    });

  smooth_corr_quotient(ngrdcol,nz,numerator,denominator,eps,corr_x_y_1);

  // Set corr_x_y_2 equal to corr_x_y_1.
  Kokkos::parallel_for("calc_comp_corrs_binormal::second_component", ngrdcol*nz,
    KOKKOS_LAMBDA(const int cell) {
      const int i = cell % ngrdcol;
      const int k = cell / ngrdcol;
      corr_x_y_2(i,k) = corr_x_y_1(i,k);
    });
}
'''
ARRAYS['smooth_corr_quotient'] = r'''
inline void smooth_corr_quotient(int ngrdcol, int nz, Real2D numerator, Real2D denominator, // In
                                 real_t denom_thresh,                                  // In
                                 Real2D quotient) {                                   // Out
  // Bound the correlation smoothly by first raising its denominator toward
  // abs(numerator)/max_mag_correlation, then toward denom_thresh. numerator
  // and denominator have the same units; quotient is dimensionless.
  Real2D tmp_denom("smooth_corr_quotient::tmp_denom",ngrdcol,nz);

  //-------------------------------------- Begin Code --------------------------------------
  Kokkos::parallel_for("smooth_corr_quotient",ngrdcol*nz,KOKKOS_LAMBDA(const int cell) {
    const int i = cell % ngrdcol;
    const int k = cell / ngrdcol;
    // First smooth_max between abs(num)/.99 and denom.
    tmp_denom(i,k) = smooth_max(Kokkos::abs(numerator(i,k))/max_mag_correlation,denominator(i,k),
                                Kokkos::min(min_max_smth_mag,denom_thresh));
    // Second smooth_max between the result and denom_thresh.
    tmp_denom(i,k) = smooth_max(tmp_denom(i,k),denom_thresh,
                                Kokkos::min(min_max_smth_mag,denom_thresh));
    quotient(i,k) = numerator(i,k)/tmp_denom(i,k);
  });
}
'''

chunks = [HEADER]
for name, info in ROUTINES.items():
    if name in ARRAYS:
        code = ARRAYS[name].replace("{DESCRIPTION}", "\n".join(description(info["block"]))).replace("{ARGUMENTS}", "\n".join(arguments(info)))
        chunks.append("\n// =============================================================================\n" + code.lstrip())
        continue
    scalar_type = "double" if name.endswith("_dp") else "real_t"
    native_head = info["block"].split("(", 1)[1].split(")", 1)[0]
    groups = [[a.strip() for a in line.replace("&", "").split(",") if a.strip()]
              for line in native_head.splitlines()]
    groups = [group for group in groups if group]
    signature = "KOKKOS_INLINE_FUNCTION " + scalar_type + " " + name + "("
    signature += ", // In\n    ".join(", ".join(scalar_type + " " + a for a in group) for group in groups) + ") {"
    lines = ["", "// =============================================================================", signature] + description(info["block"]) + ["", "  // Input Variables"] + arguments(info) + ["", "  // Return Variable", "  " + scalar_type + " " + RESULTS[name] + ";"]
    if name == "corr_LL2NN":
        lines += ["", "  // Local Variable", "  real_t log_arg; // Input into the natural logarithm [-]"]
    lines += [""] + scalar_body(name, info) + ["}"]
    chunks.append("\n".join(lines) + "\n")

# Only the three elemental whole-array calls currently used by C++ PDF callers.
for name, result_name in (("compute_variance_binormal", "xp2"), ("calc_corr_chi_x", "corr_chi_x_i"), ("calc_corr_eta_x", "corr_eta_x_i")):
    args = ROUTINES[name]["args"]
    templates = ["Field" + str(i) for i in range(len(args))]
    code = ["", "// =============================================================================", "template <" + ", ".join("class " + t for t in templates) + ">", "inline Real2D " + name + "(" + ", ".join(t + " " + a for t, a in zip(templates, args)) + ") {", "  // Apply the native elemental scalar formula to conformable execution-space", "  // column/level fields, returning a new result owned by this call.", "  static_assert(" + " && ".join(t + "::rank == 2" for t in templates) + ");", "  static_assert(" + " &&\n                ".join("Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace, typename " + t + "::memory_space>::accessible" for t in templates) + ",", '                "PDF utility arrays must be execution-space Views");', "  const int ngrdcol = " + args[0] + ".extent_int(0);", "  const int nz = " + args[0] + ".extent_int(1);", "  Real2D " + result_name + '("' + name + '::result", ngrdcol,nz);', '  Kokkos::parallel_for("' + name + '::columns",ngrdcol*nz, KOKKOS_LAMBDA(const int cell) {', "    const int i = cell % ngrdcol;", "    const int k = cell / ngrdcol;", "    " + result_name + "(i,k) = " + name + "(" + ", ".join(a + "(i,k)" for a in args) + ");", "  });", "  return " + result_name + ";", "}"]
    chunks.append("\n".join(code) + "\n")
chunks.append("\n} // namespace clubb_kokkos\n")
(OUT / "pdf_utilities.hpp").write_text("\n".join(chunks))
outline = {
    "status": "untested isolated draft",
    "canonical_source": str(SOURCE.relative_to(ROOT)),
    "canonical_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    "routines": [{"name": name, "line": info["line"], "kind": info["kind"], "arguments": info["args"], "result": RESULTS.get(name), "all_arguments_required": True} for name, info in ROUTINES.items()],
    "branch_inventory": {
        "mean_L2N": "core tiny guard; dp deliberately has none",
        "corr_NL2NN": "positive transformed stddev vs undefined fallback; ordered upper/lower clipping; dp casts core correlation limit to double",
        "corr_NN2NL": "positive transformed stddev vs original correlation; tiny guard only on ratio; ordered clipping",
        "corr_LL2NN": "both transformed stddevs positive vs original correlation; ordered clipping; core log_arg temporary, dp direct expression",
        "corr_NN2LL": "both transformed stddevs positive vs original correlation; no added denominator guard; ordered clipping",
        "calc_comp_corrs_binormal": "source numerator/denominator; smooth quotient; then second component equals first",
        "calc_corr_chi_x/calc_corr_eta_x": "strict > chi/eta tolerance; otherwise zero; then ordered clipping",
        "calc_corr_rt_x/calc_corr_thl_x": "strict > rt/thl tolerance; otherwise zero; then ordered clipping",
        "calc_xp2": "both sigma magnitudes <eps; first only; second only; neither; then negative variance->zero",
        "smooth_corr_quotient": "two sequential smooth_max calls per cell; abs(num)/max correlation then threshold; numerator divided afterward"
    },
    "dependencies": {
        "constants": "one,one_dp,zero,zero_dp,two,eps,max_mag_correlation,min_max_smth_mag,chi_tol,eta_tol,rt_tol,thl_tol; numeric_limits<real_t>::min() for Native tiny",
        "lower_native_helper": "advance_helper_module.F90:smooth_max_scalars1723–1761; exact partial draft included"
    },
    "scope_decision": "All19 public Native routines, including four dp variants and calc_xp2 omitted by JAX, are retained; no JAX gradient floors, clipping replacements or guard changes.",
    "view_adaptations": [name + " rank2 execution-space inputs -> owned Real2D" for name in ("compute_variance_binormal", "calc_corr_chi_x", "calc_corr_eta_x")],
    "bfb_specific_choices": "No new explicit FMA, forced association or math-library workaround. Native integer squares map to parenthesized products. Arithmetic equivalence/BFB pending model comparisons."
}
(OUT / "source_outline.json").write_text(json.dumps(outline, indent=2) + "\n")
print("Prepared19-routine header and Native source outline; production unchanged.")
