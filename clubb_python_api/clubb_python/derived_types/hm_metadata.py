"""Python mirror for CLUBB_core/corr_varnce_module.F90 hm_metadata_type."""

from dataclasses import dataclass, field

import numpy as np


@dataclass
class HmMetadata:
    """Python-side representation of hydrometeor metadata indices and attributes."""

    iirr: int = 0
    iirs: int = 0
    iiri: int = 0
    iirg: int = 0

    iiNr: int = 0
    iiNs: int = 0
    iiNi: int = 0
    iiNg: int = 0

    l_frozen_hm: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    l_mix_rat_hm: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    hydromet_list: list[str] = field(default_factory=list)
    hydromet_tol: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.float64))

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

    hmp2_ip_on_hmm2_ip: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.float64))
    Ncnp2_on_Ncnm2: float = 1.0
