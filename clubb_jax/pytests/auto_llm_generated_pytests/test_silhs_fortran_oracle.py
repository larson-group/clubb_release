"""Compare the inverse-normal sampler with the current Fortran routine."""

from utilities.output_paths import REPO_ROOT as ROOT
import re
import shutil
import subprocess

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.SILHS.transform_to_pdf_module import cdfnorminv
from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision

configure_jax_precision()


def test_inverse_normal_matches_source_literal_precision(tmp_path):
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("The source inverse-normal oracle requires gfortran")
    source = (ROOT / "src/SILHS/transform_to_pdf_module.F90").read_text()
    # The source routine is private. Wrap its unchanged body in a temporary
    # module, using the actual precision/constants modules; no native build or
    # fixture outputs are required by the standalone JAX model.
    routine = re.search(
        r"  subroutine cdfnorminv\(.*?  end subroutine cdfnorminv", source, re.S
    ).group()
    probabilities = np.array([
        3.0e-8, 0.001, 0.125, 0.375, 0.625, 0.875, 0.999, 1.0 - 3.0e-8,
    ])
    values = ", &\n".join(f"{value:.17e}_core_rknd" for value in probabilities)
    probe = tmp_path / "probe.f90"
    probe.write_text(
        "module native_probe\ncontains\n" + routine + "\nend module\n"
        "program probe\n"
        "use native_probe\nuse clubb_precision, only: core_rknd\n"
        "implicit none\n"
        "real(core_rknd) :: u(1,1,1,8), z(8,1,1,1)\n"
        f"u(1,1,1,:) = [ {values} ]\n"
        "call cdfnorminv(8,1,1,1,u,z)\n"
        "print '(8(ES25.17,1X))', z(:,1,1,1)\n"
        "end program\n"
    )
    executable = tmp_path / "probe"
    subprocess.run([
        compiler, "-cpp", "-DCLUBB_REAL_TYPE=8", "-O0",
        str(ROOT / "src/CLUBB_core/clubb_precision.F90"),
        str(ROOT / "src/CLUBB_core/constants_clubb.F90"),
        str(probe), "-o", str(executable),
    ], cwd=tmp_path, check=True, capture_output=True, text=True)
    reference = np.fromstring(
        subprocess.check_output([str(executable)], text=True), sep=" "
    )
    uniform = jnp.asarray(probabilities, dtype=jnp.float64).reshape(1, 1, 1, 8)
    for run in (cdfnorminv, jax.jit(cdfnorminv, static_argnums=(0, 1, 2, 3))):
        result = np.asarray(run(8, 1, 1, 1, uniform)).reshape(-1)
        # Direct double-precision literals differ by up to 8e-9. Near the
        # clipped endpoints XLA reassociation also changes log's argument;
        # check ordinary strata more tightly than those sensitive tail points.
        np.testing.assert_allclose(result[1:-1], reference[1:-1], rtol=0.0, atol=3.0e-15)
        np.testing.assert_allclose(result[[0, -1]], reference[[0, -1]], rtol=0.0, atol=2.0e-10)
