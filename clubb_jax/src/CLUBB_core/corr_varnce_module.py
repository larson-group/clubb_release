"""JAX-side port of corr_varnce_module.F90.

This module owns the default normal-space correlation arrays and the
hydrometeor PDF metadata used by hydrometeor PDF setup.

Porting deviations:
- Indices are 0-based, with -1 denoting an absent hydrometeor/PDF variable.
  Fortran uses 1-based indices and disables hydrometeors with -1 or 0,
  depending on context.
- Fortran intent(out) arguments are returned from functions.
- Fortran derived types are represented by dataclasses; ``hm_metadata_type`` is
  registered as a pytree so it can pass through JAX-transformed call graphs.
- Fortran allocation of metadata arrays is represented by NumPy arrays stored
  on immutable dataclass instances.
- Fortran debug routines that mutate ``err_info`` or write to stderr return
  booleans or formatted text in the JAX/Python port.
"""

from dataclasses import dataclass, replace

import jax
import numpy as np

from clubb_jax.src.CLUBB_core.constants_clubb import (
    eps,
    Ng_tol,
    Ni_tol,
    Nr_tol,
    Ns_tol,
    rg_tol,
    ri_tol,
    rr_tol,
    rs_tol,
)


II_CHI = 0
II_ETA = 1
II_W = 2
II_NCN = 3
II_RR = 4
II_NR = 5
II_RI = 6
II_NI = 7
II_RS = 8
II_NS = 9
II_RG = 10
II_NG = 11


@dataclass
class HmMetadata:
    """Small compatibility metadata object used by unit tests and index helpers."""

    hydromet_dim: int = 0
    iirr: int = -1
    iirs: int = -1
    iiri: int = -1
    iirg: int = -1
    iiNr: int = -1
    iiNs: int = -1
    iiNi: int = -1
    iiNg: int = -1
    iiPDF_chi: int = II_CHI
    iiPDF_eta: int = II_ETA
    iiPDF_w: int = II_W
    iiPDF_Ncn: int = II_NCN
    iiPDF_rr: int = -1
    iiPDF_rs: int = -1
    iiPDF_ri: int = -1
    iiPDF_rg: int = -1
    iiPDF_Nr: int = -1
    iiPDF_Ns: int = -1
    iiPDF_Ni: int = -1
    iiPDF_Ng: int = -1

    @property
    def pdf_dim(self):
        return max(
            self.iiPDF_chi,
            self.iiPDF_eta,
            self.iiPDF_w,
            self.iiPDF_Ncn,
            self.iiPDF_rr,
            self.iiPDF_rs,
            self.iiPDF_ri,
            self.iiPDF_rg,
            self.iiPDF_Nr,
            self.iiPDF_Ns,
            self.iiPDF_Ni,
            self.iiPDF_Ng,
        ) + 1


@dataclass(frozen=True)
class hmp2_ip_on_hmm2_ip_ratios_type:
    """Prescribed parameters for hydrometeor values of <hm|_ip'^2> / <hm|_ip>^2,
    where  <hm|_ip'^2>  is the in-precip. variance of the hydrometeor and
           <hm|_ip>     is the in-precip. mean of the hydrometeor

    These values are dependent on the horizontal grid spacing of the run, and are calculated
    using a slope and intercept corresponding to each hydrometer.
    """
    rr: float = 1.0
    Nr: float = 1.0
    ri: float = 1.0
    Ni: float = 1.0
    rs: float = 1.0
    Ns: float = 1.0
    rg: float = 1.0
    Ng: float = 1.0


@dataclass(frozen=True)
class hmp2_ip_on_hmm2_ip_slope_type:
    """These slopes and intercepts below are used to calculate the hmp2_ip_on_hmm2_ip_ratios_type
    values that are defined above. This functionality is described by equations 8, 10, & 11
    in ``Parameterization of the Spatial Variability of Rain for Large-Scale Models and
    Remote Sensing`` Lebo, et al, October 2015
    https://journals.ametsoc.org/doi/pdf/10.1175/JAMC-D-15-0066.1
    see clubb:ticket:830 for more detail

    hmp2_ip_on_hmm2_ip(iirr) = hmp2_ip_on_hmm2_ip_intrcpt%rr +
                               hmp2_ip_on_hmm2_ip_slope%rr * max( host_dx, host_dy )

    In Lebo et al. the suggested values were
        slope = 2.12e-5 [1/m]
        intercept = 0.54 [-]

    In CLUBB standalone, these parameters can be set based on the value for a
    given case in the CASE_model.in file.
    """
    rr: float = 2.12e-5
    Nr: float = 2.12e-5
    ri: float = 2.12e-5
    Ni: float = 2.12e-5
    rs: float = 2.12e-5
    Ns: float = 2.12e-5
    rg: float = 2.12e-5
    Ng: float = 2.12e-5


@dataclass(frozen=True)
class hmp2_ip_on_hmm2_ip_intrcpt_type:
    rr: float = 0.54
    Nr: float = 0.54
    ri: float = 0.54
    Ni: float = 0.54
    rs: float = 0.54
    Ns: float = 0.54
    rg: float = 0.54
    Ng: float = 0.54


@jax.tree_util.register_pytree_node_class
@dataclass(frozen=True)
class hm_metadata_type:
    """Hydrometeor metadata.

    Variables
    Microphysics mixing ratios
    Microphysics concentrations
    Logical fields
    Latin hypercube indices / Correlation array indices
    """
    iirr: int = -1
    iirs: int = -1
    iiri: int = -1
    iirg: int = -1
    iiNr: int = -1
    iiNs: int = -1
    iiNi: int = -1
    iiNg: int = -1
    l_frozen_hm: np.ndarray | None = None
    l_mix_rat_hm: np.ndarray | None = None
    hydromet_list: tuple[str, ...] = ()
    hydromet_tol: np.ndarray | None = None
    iiPDF_chi: int = -1
    iiPDF_eta: int = -1
    iiPDF_w: int = -1
    iiPDF_rr: int = -1
    iiPDF_rs: int = -1
    iiPDF_ri: int = -1
    iiPDF_rg: int = -1
    iiPDF_Nr: int = -1
    iiPDF_Ns: int = -1
    iiPDF_Ni: int = -1
    iiPDF_Ng: int = -1
    iiPDF_Ncn: int = -1
    hmp2_ip_on_hmm2_ip: np.ndarray | None = None
    Ncnp2_on_Ncnm2: float = 1.0

    @property
    def hydromet_dim(self):
        if self.hydromet_list:
            return len(self.hydromet_list)
        if self.hmp2_ip_on_hmm2_ip is not None:
            return len(self.hmp2_ip_on_hmm2_ip)
        return max(
            self.iirr,
            self.iirs,
            self.iiri,
            self.iirg,
            self.iiNr,
            self.iiNs,
            self.iiNi,
            self.iiNg,
        ) + 1

    @property
    def pdf_dim(self):
        return max(
            self.iiPDF_chi,
            self.iiPDF_eta,
            self.iiPDF_w,
            self.iiPDF_rr,
            self.iiPDF_rs,
            self.iiPDF_ri,
            self.iiPDF_rg,
            self.iiPDF_Nr,
            self.iiPDF_Ns,
            self.iiPDF_Ni,
            self.iiPDF_Ng,
            self.iiPDF_Ncn,
        ) + 1

    def tree_flatten(self):
        children = (
            self.l_frozen_hm,
            self.l_mix_rat_hm,
            self.hydromet_tol,
            self.hmp2_ip_on_hmm2_ip,
        )
        aux_data = (
            self.iirr,
            self.iirs,
            self.iiri,
            self.iirg,
            self.iiNr,
            self.iiNs,
            self.iiNi,
            self.iiNg,
            self.hydromet_list,
            self.iiPDF_chi,
            self.iiPDF_eta,
            self.iiPDF_w,
            self.iiPDF_rr,
            self.iiPDF_rs,
            self.iiPDF_ri,
            self.iiPDF_rg,
            self.iiPDF_Nr,
            self.iiPDF_Ns,
            self.iiPDF_Ni,
            self.iiPDF_Ng,
            self.iiPDF_Ncn,
            self.Ncnp2_on_Ncnm2,
        )
        return children, aux_data

    @classmethod
    def tree_unflatten(cls, aux_data, children):
        (
            iirr,
            iirs,
            iiri,
            iirg,
            iiNr,
            iiNs,
            iiNi,
            iiNg,
            hydromet_list,
            iiPDF_chi,
            iiPDF_eta,
            iiPDF_w,
            iiPDF_rr,
            iiPDF_rs,
            iiPDF_ri,
            iiPDF_rg,
            iiPDF_Nr,
            iiPDF_Ns,
            iiPDF_Ni,
            iiPDF_Ng,
            iiPDF_Ncn,
            Ncnp2_on_Ncnm2,
        ) = aux_data
        (
            l_frozen_hm,
            l_mix_rat_hm,
            hydromet_tol,
            hmp2_ip_on_hmm2_ip,
        ) = children
        return cls(
            iirr=iirr,
            iirs=iirs,
            iiri=iiri,
            iirg=iirg,
            iiNr=iiNr,
            iiNs=iiNs,
            iiNi=iiNi,
            iiNg=iiNg,
            l_frozen_hm=l_frozen_hm,
            l_mix_rat_hm=l_mix_rat_hm,
            hydromet_list=hydromet_list,
            hydromet_tol=hydromet_tol,
            iiPDF_chi=iiPDF_chi,
            iiPDF_eta=iiPDF_eta,
            iiPDF_w=iiPDF_w,
            iiPDF_rr=iiPDF_rr,
            iiPDF_rs=iiPDF_rs,
            iiPDF_ri=iiPDF_ri,
            iiPDF_rg=iiPDF_rg,
            iiPDF_Nr=iiPDF_Nr,
            iiPDF_Ns=iiPDF_Ns,
            iiPDF_Ni=iiPDF_Ni,
            iiPDF_Ng=iiPDF_Ng,
            iiPDF_Ncn=iiPDF_Ncn,
            hmp2_ip_on_hmm2_ip=hmp2_ip_on_hmm2_ip,
            Ncnp2_on_Ncnm2=Ncnp2_on_Ncnm2,
        )


d_var_total = 12

# Size of the default correlation arrays
corr_array_n_cloud_def = np.array(
    [
        #  chi   eta    w      Ncn     rr      Nr      ri      Ni      rs      Ns      rg      Ng
        1.0, -0.6, 0.09, 0.09, 0.788, 0.675, 0.240, 0.222, 0.240, 0.222, 0.240, 0.222,
        0.0, 1.0, 0.027, 0.027, 0.114, 0.115, -0.029, 0.093, 0.022, 0.013, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.34, 0.315, 0.270, 0.120, 0.167, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.464, 0.320, 0.168, 0.232, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 1.0, 0.821, 0.0, 0.0, 0.173, 0.164, 0.319, 0.308,
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.152, 0.143, 0.0, 0.0, 0.285, 0.273,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.758, 0.585, 0.571, 0.379, 0.363,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.571, 0.550, 0.363, 0.345,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.758, 0.485, 0.470,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.470, 0.450,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.758,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ],
    dtype=np.float64,
).reshape((d_var_total, d_var_total), order="F")

corr_array_n_below_def = np.array(
    [
        #  chi   eta    w      Ncn     rr      Nr      ri      Ni      rs      Ns      rg      Ng
        1.0, 0.3, 0.09, 0.09, 0.788, 0.675, 0.240, 0.222, 0.240, 0.222, 0.240, 0.222,
        0.0, 1.0, 0.027, 0.027, 0.114, 0.115, -0.029, 0.093, 0.022, 0.013, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.34, 0.315, 0.270, 0.120, 0.167, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.464, 0.320, 0.168, 0.232, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 1.0, 0.821, 0.0, 0.0, 0.173, 0.164, 0.319, 0.308,
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.152, 0.143, 0.0, 0.0, 0.285, 0.273,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.758, 0.585, 0.571, 0.379, 0.363,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.571, 0.550, 0.363, 0.345,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.758, 0.485, 0.470,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.470, 0.450,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.758,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ],
    dtype=np.float64,
).reshape((d_var_total, d_var_total), order="F")

CORR_ARRAY_N_CLOUD_DEF = corr_array_n_cloud_def
CORR_ARRAY_N_BELOW_DEF = corr_array_n_below_def


def def_corr_idx(iiPDF_x, hm_metadata):
    """Map from a iiPDF index to the corresponding index in the default
    correlation arrays."""
    # Constant Parameters

    # Indices that represent the order in the default corr arrays
    # (chi (old s), eta (old t), w, Ncn, rr, Nr, ri, Ni, rs, Ns, rg, Ng)
    ii_chi = 0
    ii_eta = 1
    ii_w = 2
    ii_Ncn = 3
    ii_rr = 4
    ii_Nr = 5
    ii_ri = 6
    ii_Ni = 7
    ii_rs = 8
    ii_Ns = 9
    ii_rg = 10
    ii_Ng = 11

    ii_def_corr = -1

    if iiPDF_x == hm_metadata.iiPDF_chi:
        ii_def_corr = ii_chi
    elif iiPDF_x == hm_metadata.iiPDF_eta:
        ii_def_corr = ii_eta
    elif iiPDF_x == hm_metadata.iiPDF_w:
        ii_def_corr = ii_w
    elif iiPDF_x == hm_metadata.iiPDF_Ncn:
        ii_def_corr = ii_Ncn
    elif iiPDF_x == hm_metadata.iiPDF_rr:
        ii_def_corr = ii_rr
    elif iiPDF_x == hm_metadata.iiPDF_Nr:
        ii_def_corr = ii_Nr
    elif iiPDF_x == hm_metadata.iiPDF_ri:
        ii_def_corr = ii_ri
    elif iiPDF_x == hm_metadata.iiPDF_Ni:
        ii_def_corr = ii_Ni
    elif iiPDF_x == hm_metadata.iiPDF_rs:
        ii_def_corr = ii_rs
    elif iiPDF_x == hm_metadata.iiPDF_Ns:
        ii_def_corr = ii_Ns
    elif iiPDF_x == hm_metadata.iiPDF_rg:
        ii_def_corr = ii_rg
    elif iiPDF_x == hm_metadata.iiPDF_Ng:
        ii_def_corr = ii_Ng

    return ii_def_corr


def _def_corr_idx_from(iiPDF_x, hm_metadata):
    if isinstance(hm_metadata, (tuple, list, np.ndarray)):
        return int(hm_metadata[iiPDF_x])
    return def_corr_idx(iiPDF_x, hm_metadata)


def set_corr_arrays_to_default(pdf_dim, hm_metadata):
    """If there are no corr_array.in files for the current case, default
    correlations are used."""
    corr_array_n_cloud = np.zeros((pdf_dim, pdf_dim), dtype=np.float64)
    corr_array_n_below = np.zeros((pdf_dim, pdf_dim), dtype=np.float64)

    for i in range(pdf_dim):
        corr_array_n_cloud[i, i] = 1.0
        corr_array_n_below[i, i] = 1.0

    for i in range(pdf_dim - 1):
        for j in range(i + 1, pdf_dim):
            idx_i = _def_corr_idx_from(i, hm_metadata)
            idx_j = _def_corr_idx_from(j, hm_metadata)

            if idx_i > idx_j:
                corr_array_n_cloud[j, i] = corr_array_n_cloud_def[idx_i, idx_j]
                corr_array_n_below[j, i] = corr_array_n_below_def[idx_i, idx_j]
            else:
                corr_array_n_cloud[j, i] = corr_array_n_cloud_def[idx_j, idx_i]
                corr_array_n_below[j, i] = corr_array_n_below_def[idx_j, idx_i]

    return corr_array_n_cloud, corr_array_n_below


def kk_prescribed_correlations():
    """Return prescribed normal-space correlations for the warm-rain KK PDF."""
    return {
        "corr_chi_eta": corr_array_n_cloud_def[II_ETA, II_CHI],
        "corr_w_chi": corr_array_n_cloud_def[II_W, II_CHI],
        "corr_w_eta": corr_array_n_cloud_def[II_W, II_ETA],
        "corr_w_Ncn": corr_array_n_cloud_def[II_NCN, II_W],
        "corr_w_rr": corr_array_n_cloud_def[II_RR, II_W],
        "corr_w_Nr": corr_array_n_cloud_def[II_NR, II_W],
        "corr_chi_Ncn": corr_array_n_cloud_def[II_NCN, II_CHI],
        "corr_chi_rr": corr_array_n_cloud_def[II_RR, II_CHI],
        "corr_chi_Nr": corr_array_n_cloud_def[II_NR, II_CHI],
        "corr_eta_Ncn": corr_array_n_cloud_def[II_NCN, II_ETA],
        "corr_eta_rr": corr_array_n_cloud_def[II_RR, II_ETA],
        "corr_eta_Nr": corr_array_n_cloud_def[II_NR, II_ETA],
        "corr_rr_Nr": corr_array_n_cloud_def[II_NR, II_RR],
    }


def get_corr_var_index(var_name, hm_metadata):
    """Returns the index for a variable based on its name."""
    i = -1

    var_name = var_name.strip()
    if var_name == "chi":
        i = hm_metadata.iiPDF_chi
    elif var_name == "eta":
        i = hm_metadata.iiPDF_eta
    elif var_name == "w":
        i = hm_metadata.iiPDF_w
    elif var_name == "Ncn":
        i = hm_metadata.iiPDF_Ncn
    elif var_name == "rr":
        i = hm_metadata.iiPDF_rr
    elif var_name == "Nr":
        i = hm_metadata.iiPDF_Nr
    elif var_name == "ri":
        i = hm_metadata.iiPDF_ri
    elif var_name == "Ni":
        i = hm_metadata.iiPDF_Ni
    elif var_name == "rs":
        i = hm_metadata.iiPDF_rs
    elif var_name == "Ns":
        i = hm_metadata.iiPDF_Ns
    elif var_name == "rg":
        i = hm_metadata.iiPDF_rg
    elif var_name == "Ng":
        i = hm_metadata.iiPDF_Ng

    return i


def init_pdf_hydromet_arrays_api(
    host_dx,
    host_dy,
    hydromet_dim,
    iirr,
    iiNr,
    iiri,
    iiNi,
    iirs,
    iiNs,
    iirg,
    iiNg,
    Ncnp2_on_Ncnm2,
    hmp2_ip_on_hmm2_ip_slope_in=None,
    hmp2_ip_on_hmm2_ip_intrcpt_in=None,
):
    """Initialize hydrometeor metadata and return ``hm_metadata, pdf_dim``.

    DESCRIPTION:
        This subroutine intializes the hydromet arrays(iirr, iiNr, etc.) to the values specified by
        the input arguements, this determines which hyrometeors are to be used by the microphysics
        scheme. It also sets up the corresponding pdf and hydromet arrays, and calculates the
        subgrid variance ratio for each hydrometeor.

    OPTIONAL FUNCTIONALITY:
        The subgrid variance ratio for each hydrometeor is calculated based on the grid spacing
        defined by the host model. The calculation is a linear equation defined by a slope and
        intercept, each of which may or may not be passed in to this subroutine. If the slope
        and/or intercept are not passed in through the arguement list the default values, which
        are set in the corresponding type definitions, will be used. Otherwise the values
        specified by the aruements will be used.

    NOTES:
        'hmp2_ip_on_hmm2_ip_slope_in' is of type 'hmp2_ip_on_hmm2_ip_slope_type' and
        'hmp2_ip_on_hmm2_ip_intrcpt_in' is of type 'hmp2_ip_on_hmm2_ip_intrcpt_in', both of which
        are deinfed in corr_vrnce_module.F90, and made public through this API.

        If full control over the hydrometeor variance ratios is desired, pass in slopes that are
        initialized to 0.0, this causes the ratios to no longer depend on the grid spacing. Then
        pass in the intercepts set to the values of the desired ratios.
    """
    hm_metadata = hm_metadata_type(
        iirr=iirr,
        iirs=iirs,
        iiri=iiri,
        iirg=iirg,
        iiNr=iiNr,
        iiNs=iiNs,
        iiNi=iiNi,
        iiNg=iiNg,
        Ncnp2_on_Ncnm2=Ncnp2_on_Ncnm2,
    )

    hmp2_ip_on_hmm2_ip_slope = (
        hmp2_ip_on_hmm2_ip_slope_in
        if hmp2_ip_on_hmm2_ip_slope_in is not None
        else hmp2_ip_on_hmm2_ip_slope_type()
    )
    hmp2_ip_on_hmm2_ip_intrcpt = (
        hmp2_ip_on_hmm2_ip_intrcpt_in
        if hmp2_ip_on_hmm2_ip_intrcpt_in is not None
        else hmp2_ip_on_hmm2_ip_intrcpt_type()
    )

    #-----------------------------------------------------------
    # Calculate the subgrid variances of the hydrometeors
    #-----------------------------------------------------------

    # If slope and intercept are present in call, then overwrite default values
    max_host_delta = max(host_dx, host_dy)
    hmp2_ip_on_hmm2_ip = np.zeros(hydromet_dim, dtype=np.float64)

    if iirr >= 0:
        hmp2_ip_on_hmm2_ip[iirr] = (
            hmp2_ip_on_hmm2_ip_intrcpt.rr
            + hmp2_ip_on_hmm2_ip_slope.rr * max_host_delta
        )
    if iirs >= 0:
        hmp2_ip_on_hmm2_ip[iirs] = (
            hmp2_ip_on_hmm2_ip_intrcpt.rs
            + hmp2_ip_on_hmm2_ip_slope.rs * max_host_delta
        )
    if iiri >= 0:
        hmp2_ip_on_hmm2_ip[iiri] = (
            hmp2_ip_on_hmm2_ip_intrcpt.ri
            + hmp2_ip_on_hmm2_ip_slope.ri * max_host_delta
        )
    if iirg >= 0:
        hmp2_ip_on_hmm2_ip[iirg] = (
            hmp2_ip_on_hmm2_ip_intrcpt.rg
            + hmp2_ip_on_hmm2_ip_slope.rg * max_host_delta
        )
    if iiNr >= 0:
        hmp2_ip_on_hmm2_ip[iiNr] = (
            hmp2_ip_on_hmm2_ip_intrcpt.Nr
            + hmp2_ip_on_hmm2_ip_slope.Nr * max_host_delta
        )
    if iiNs >= 0:
        hmp2_ip_on_hmm2_ip[iiNs] = (
            hmp2_ip_on_hmm2_ip_intrcpt.Ns
            + hmp2_ip_on_hmm2_ip_slope.Ns * max_host_delta
        )
    if iiNi >= 0:
        hmp2_ip_on_hmm2_ip[iiNi] = (
            hmp2_ip_on_hmm2_ip_intrcpt.Ni
            + hmp2_ip_on_hmm2_ip_slope.Ni * max_host_delta
        )
    if iiNg >= 0:
        hmp2_ip_on_hmm2_ip[iiNg] = (
            hmp2_ip_on_hmm2_ip_intrcpt.Ng
            + hmp2_ip_on_hmm2_ip_slope.Ng * max_host_delta
        )

    hydromet_list = [""] * hydromet_dim
    hydromet_tol = np.zeros(hydromet_dim, dtype=np.float64)
    l_mix_rat_hm = np.zeros(hydromet_dim, dtype=bool)
    l_frozen_hm = np.zeros(hydromet_dim, dtype=bool)

    if iirr >= 0:
        # The microphysics scheme predicts rain water mixing ratio, rr.
        hydromet_list[iirr] = "rrm"
        l_mix_rat_hm[iirr] = True
        l_frozen_hm[iirr] = False
        hydromet_tol[iirr] = rr_tol
    if iiri >= 0:
        # The microphysics scheme predicts ice mixing ratio, ri.
        hydromet_list[iiri] = "rim"
        l_mix_rat_hm[iiri] = True
        l_frozen_hm[iiri] = True
        hydromet_tol[iiri] = ri_tol
    if iirs >= 0:
        # The microphysics scheme predicts snow mixing ratio, rs.
        hydromet_list[iirs] = "rsm"
        l_mix_rat_hm[iirs] = True
        l_frozen_hm[iirs] = True
        hydromet_tol[iirs] = rs_tol
    if iirg >= 0:
        # The microphysics scheme predicts graupel mixing ratio, rg.
        hydromet_list[iirg] = "rgm"
        l_mix_rat_hm[iirg] = True
        l_frozen_hm[iirg] = True
        hydromet_tol[iirg] = rg_tol
    if iiNr >= 0:
        # The microphysics scheme predicts rain drop concentration, Nr.
        hydromet_list[iiNr] = "Nrm"
        l_frozen_hm[iiNr] = False
        l_mix_rat_hm[iiNr] = False
        hydromet_tol[iiNr] = Nr_tol
    if iiNi >= 0:
        # The microphysics scheme predicts ice concentration, Ni.
        hydromet_list[iiNi] = "Nim"
        l_mix_rat_hm[iiNi] = False
        l_frozen_hm[iiNi] = True
        hydromet_tol[iiNi] = Ni_tol
    if iiNs >= 0:
        # The microphysics scheme predicts snowflake concentration, Ns.
        hydromet_list[iiNs] = "Nsm"
        l_mix_rat_hm[iiNs] = False
        l_frozen_hm[iiNs] = True
        hydromet_tol[iiNs] = Ns_tol
    if iiNg >= 0:
        # The microphysics scheme predicts graupel concentration, Ng.
        hydromet_list[iiNg] = "Ngm"
        l_mix_rat_hm[iiNg] = False
        l_frozen_hm[iiNg] = True
        hydromet_tol[iiNg] = Ng_tol

    hm_metadata = replace(
        hm_metadata,
        hmp2_ip_on_hmm2_ip=hmp2_ip_on_hmm2_ip,
        hydromet_list=tuple(hydromet_list),
        hydromet_tol=hydromet_tol,
        l_mix_rat_hm=l_mix_rat_hm,
        l_frozen_hm=l_frozen_hm,
        iiPDF_chi=0,
        iiPDF_eta=1,
        iiPDF_w=2,
        iiPDF_Ncn=3,
    )

    #-----------------------------------------------------------
    # Set up the PDF indices
    #-----------------------------------------------------------

    pdf_count = hm_metadata.iiPDF_Ncn

    if hydromet_dim > 0:
        # Loop over hydrometeors.
        # Hydrometeor indices in the PDF arrays should be in the same order as
        # found in the hydrometeor arrays.
        for i in range(hydromet_dim):
            if i == iirr:
                pdf_count += 1
                hm_metadata = replace(hm_metadata, iiPDF_rr=pdf_count)
            if i == iiNr:
                pdf_count += 1
                hm_metadata = replace(hm_metadata, iiPDF_Nr=pdf_count)
            if i == iiri:
                pdf_count += 1
                hm_metadata = replace(hm_metadata, iiPDF_ri=pdf_count)
            if i == iiNi:
                pdf_count += 1
                hm_metadata = replace(hm_metadata, iiPDF_Ni=pdf_count)
            if i == iirs:
                pdf_count += 1
                hm_metadata = replace(hm_metadata, iiPDF_rs=pdf_count)
            if i == iiNs:
                pdf_count += 1
                hm_metadata = replace(hm_metadata, iiPDF_Ns=pdf_count)
            if i == iirg:
                pdf_count += 1
                hm_metadata = replace(hm_metadata, iiPDF_rg=pdf_count)
            if i == iiNg:
                pdf_count += 1
                hm_metadata = replace(hm_metadata, iiPDF_Ng=pdf_count)

    pdf_dim = pdf_count + 1

    return hm_metadata, pdf_dim


def init_pdf_hydromet_arrays(
    host_dx=0.0,
    host_dy=0.0,
    hydromet_dim=0,
    iirr=-1,
    iiNr=-1,
    iiri=-1,
    iiNi=-1,
    iirs=-1,
    iiNs=-1,
    iirg=-1,
    iiNg=-1,
    Ncnp2_on_Ncnm2=1.0,
):
    """Compatibility wrapper returning only hydrometeor metadata."""
    hm_metadata, _ = init_pdf_hydromet_arrays_api(
        host_dx,
        host_dy,
        hydromet_dim,
        iirr,
        iiNr,
        iiri,
        iiNi,
        iirs,
        iiNs,
        iirg,
        iiNg,
        Ncnp2_on_Ncnm2,
    )
    return hm_metadata


def kk_hm_metadata():
    """Return the warm-rain KK hydrometeor metadata used by the RICO-style setup."""
    hm_metadata = init_pdf_hydromet_arrays(
        hydromet_dim=2,
        iirr=0,
        iiNr=1,
        host_dx=0.0,
        host_dy=0.0,
    )
    return replace(
        hm_metadata,
        hmp2_ip_on_hmm2_ip=np.full(2, 1.25, dtype=np.float64),
    )


def setup_corr_varnce_array_api(
    pdf_dim,
    hm_metadata,
    l_fix_w_chi_eta_correlations,
    corr_array_n_cloud_default=None,
    corr_array_n_below_default=None,
):
    """Set up full symmetric cloud and below-cloud normal-space correlations."""
    del l_fix_w_chi_eta_correlations

    if corr_array_n_cloud_default is not None and corr_array_n_below_default is not None:
        corr_array_n_cloud = np.array(corr_array_n_cloud_default, dtype=np.float64)
        corr_array_n_below = np.array(corr_array_n_below_default, dtype=np.float64)
    else:
        corr_array_n_cloud, corr_array_n_below = set_corr_arrays_to_default(
            pdf_dim,
            hm_metadata,
        )

    corr_array_n_cloud = (
        np.tril(corr_array_n_cloud)
        + np.tril(corr_array_n_cloud, -1).T
    )
    corr_array_n_below = (
        np.tril(corr_array_n_below)
        + np.tril(corr_array_n_below, -1).T
    )

    return corr_array_n_cloud, corr_array_n_below


def assert_corr_symmetric(corr_array_n):
    """Return True when a normal-space correlation matrix is symmetric with unit diagonal."""
    tol = 1.0e-6
    corr_array_n = np.asarray(corr_array_n, dtype=np.float64)

    symmetric = np.all(np.abs(corr_array_n - corr_array_n.T) <= tol)
    unit_diagonal = np.all(np.abs(np.diagonal(corr_array_n) - 1.0) <= eps)

    return bool(symmetric and unit_diagonal)


def print_corr_matrix(pdf_dim, corr_array_n, stream=None):
    """Print a correlation matrix in the same orientation as the Fortran debug routine."""
    output = []

    for n in range(pdf_dim):
        row = []
        for m in range(pdf_dim):
            row.append(f"{corr_array_n[m, n]:5.2f}")
        output.append(" ".join(row))

    text = "\n".join(output)
    if stream is None:
        print(text)
    else:
        print(text, file=stream)
