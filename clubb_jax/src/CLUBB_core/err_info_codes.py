"""JAX-side error reason codes for `ErrInfo`."""

from __future__ import annotations

import numpy as np

ERR_NONE = 0
ERR_XP2_XPYP_MULTIPLE_LHS_REQUIRED = 1
ERR_XP2_XPYP_INVALID_C_UU = 2
ERR_RADIATION_COS_SOLAR_ZEN_TIME_NOT_FOUND = 3
ERR_RADIATION_COS_SOLAR_ZEN_INVALID_HOUR = 4


ERROR_MESSAGES = {
    ERR_NONE: "",
    ERR_XP2_XPYP_MULTIPLE_LHS_REQUIRED: (
        "advance_xp2_xpyp requires the multiple-LHS solve because C2rt, "
        "C2thl, and C2rtthl are not equivalent."
    ),
    ERR_XP2_XPYP_INVALID_C_UU: (
        "advance_xp2_xpyp found invalid C_uu_shr or C_uu_buoy values."
    ),
    ERR_RADIATION_COS_SOLAR_ZEN_TIME_NOT_FOUND: (
        "radiation_driver did not find time_current in cos_solar_zen_times."
    ),
    ERR_RADIATION_COS_SOLAR_ZEN_INVALID_HOUR: (
        "cos_solar_zen received an invalid current_time."
    ),
}


def message_for(code: int) -> str:
    code = int(code)
    return ERROR_MESSAGES.get(code, f"Unknown CLUBB JAX error reason code: {code}")


def nonzero_codes(reason_code) -> list[int]:
    if reason_code is None:
        return []
    codes = np.asarray(reason_code, dtype=np.int32).reshape(-1)
    return sorted({int(code) for code in codes if int(code) != ERR_NONE})


def messages_for(reason_code) -> list[str]:
    return [message_for(code) for code in nonzero_codes(reason_code)]
