"""Shared probability/cloud-transport coloring for PDF teaching plots."""

from __future__ import annotations

from functools import lru_cache

import numpy as np


BIVARIATE_LEVELS = 16


def robust_field_upper(fields) -> float:
    """Return a shared robust positive upper scale for one or more fields."""
    positive = []
    for field in fields:
        values = np.nan_to_num(
            np.asarray(field, dtype=float), nan=0.0, posinf=0.0, neginf=0.0
        )
        values = values[values > 0.0]
        if values.size:
            positive.append(values)
    if not positive:
        return 0.0
    return float(np.percentile(np.concatenate(positive), 99.5))


def normalize_reference_field(field: np.ndarray, upper: float | None = None) -> np.ndarray:
    """Log-normalize a positive field without letting one sharp peak dominate."""
    field = np.nan_to_num(
        np.asarray(field, dtype=float), nan=0.0, posinf=0.0, neginf=0.0
    )
    if upper is None:
        upper = robust_field_upper((field,))
    upper = max(float(upper), 0.0)
    if upper <= 0.0:
        return np.zeros_like(field)
    knee = max(0.025 * upper, np.finfo(float).tiny)
    normalized = np.log1p(np.maximum(field, 0.0) / knee) / np.log1p(upper / knee)
    return np.clip(normalized, 0.0, 1.0)


def bivariate_reference_rgb(
    probability: float,
    transport: float,
    transport_color=(37.0, 99.0, 235.0),
) -> tuple[int, int, int]:
    """Map probability to gold, a transport signal to color, and overlap to pale."""
    probability = float(np.clip(probability, 0.0, 1.0))
    transport = float(np.clip(transport, 0.0, 1.0))
    dark = np.array([15.0, 23.0, 42.0])
    probability_color = np.array([245.0, 158.0, 11.0])
    transport_color = np.asarray(transport_color, dtype=float)
    overlap_color = np.array([255.0, 247.0, 204.0])
    strength = 1.0 - (1.0 - probability) * (1.0 - transport)
    mixture = transport / max(probability + transport, 1.0e-12)
    hue = (1.0 - mixture) * probability_color + mixture * transport_color
    color = (1.0 - strength) * dark + strength * hue
    overlap = min(probability, transport) ** 0.75
    color = (1.0 - 0.72 * overlap) * color + 0.72 * overlap * overlap_color
    return tuple(int(round(value)) for value in np.clip(color, 0.0, 255.0))


@lru_cache(maxsize=16)
def bivariate_reference_colorscale(
    levels: int = BIVARIATE_LEVELS,
    transport_color=(37.0, 99.0, 235.0),
):
    """Return the quantized Plotly colorscale paired with ``bivariate_codes``."""
    colorscale = []
    count = levels * levels
    for probability_bin in range(levels):
        for transport_bin in range(levels):
            index = probability_bin * levels + transport_bin
            probability = probability_bin / (levels - 1)
            transport = transport_bin / (levels - 1)
            red, green, blue = bivariate_reference_rgb(
                probability,
                transport,
                transport_color,
            )
            colorscale.append([index / (count - 1), f"rgb({red},{green},{blue})"])
    return colorscale


def bivariate_codes(
    probability,
    transport,
    probability_upper: float | None = None,
    transport_upper: float | None = None,
    levels: int = BIVARIATE_LEVELS,
):
    """Encode two positive fields and return codes plus their normalized values."""
    normalized_probability = normalize_reference_field(probability, probability_upper)
    normalized_transport = normalize_reference_field(transport, transport_upper)
    probability_bin = np.rint(normalized_probability * (levels - 1)).astype(int)
    transport_bin = np.rint(normalized_transport * (levels - 1)).astype(int)
    return (
        probability_bin * levels + transport_bin,
        normalized_probability,
        normalized_transport,
    )


def signed_bivariate_reference_rgb(
    probability: float,
    signed_transport: float,
) -> tuple[int, int, int]:
    """Map probability to gold and signed cloud transport to red/blue.

    Positive transport is blue, negative transport is red, and probability is
    gold.  Where probability and transport are both strong, the color is
    lightened without discarding the transport sign.
    """
    probability = float(np.clip(probability, 0.0, 1.0))
    signed_transport = float(np.clip(signed_transport, -1.0, 1.0))
    transport = abs(signed_transport)
    dark = np.array([15.0, 23.0, 42.0])
    probability_color = np.array([245.0, 158.0, 11.0])
    if signed_transport >= 0.0:
        transport_color = np.array([37.0, 99.0, 235.0])
        overlap_color = np.array([219.0, 234.0, 254.0])
    else:
        transport_color = np.array([239.0, 68.0, 68.0])
        overlap_color = np.array([254.0, 205.0, 211.0])
    strength = 1.0 - (1.0 - probability) * (1.0 - transport)
    mixture = transport / max(probability + transport, 1.0e-12)
    hue = (1.0 - mixture) * probability_color + mixture * transport_color
    color = (1.0 - strength) * dark + strength * hue
    overlap = min(probability, transport) ** 0.75
    color = (1.0 - 0.68 * overlap) * color + 0.68 * overlap * overlap_color
    return tuple(int(round(value)) for value in np.clip(color, 0.0, 255.0))


@lru_cache(maxsize=8)
def signed_bivariate_reference_colorscale(levels: int = BIVARIATE_LEVELS):
    """Return the quantized colorscale paired with ``signed_bivariate_codes``."""
    signed_levels = 2 * levels - 1
    count = levels * signed_levels
    colorscale = []
    for probability_bin in range(levels):
        for signed_bin in range(signed_levels):
            index = probability_bin * signed_levels + signed_bin
            probability = probability_bin / (levels - 1)
            signed_transport = (signed_bin - (levels - 1)) / (levels - 1)
            red, green, blue = signed_bivariate_reference_rgb(
                probability, signed_transport
            )
            colorscale.append([index / (count - 1), f"rgb({red},{green},{blue})"])
    return colorscale


def signed_bivariate_codes(
    probability,
    signed_transport,
    probability_upper: float | None = None,
    transport_upper: float | None = None,
    levels: int = BIVARIATE_LEVELS,
):
    """Encode probability and a signed transport field into one color index."""
    probability = np.asarray(probability, dtype=float)
    signed_transport = np.nan_to_num(
        np.asarray(signed_transport, dtype=float), nan=0.0, posinf=0.0, neginf=0.0
    )
    normalized_probability = normalize_reference_field(probability, probability_upper)
    normalized_magnitude = normalize_reference_field(
        np.abs(signed_transport), transport_upper
    )
    normalized_signed = np.sign(signed_transport) * normalized_magnitude
    probability_bin = np.rint(normalized_probability * (levels - 1)).astype(int)
    signed_bin = np.rint(normalized_signed * (levels - 1)).astype(int) + levels - 1
    signed_levels = 2 * levels - 1
    return (
        probability_bin * signed_levels + signed_bin,
        normalized_probability,
        normalized_signed,
    )
