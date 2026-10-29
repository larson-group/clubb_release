import numpy as np

from dash_app.shared.bivariate_heatmap import (
    BIVARIATE_LEVELS,
    signed_bivariate_codes,
    signed_bivariate_reference_colorscale,
    signed_bivariate_reference_rgb,
)


def test_signed_bivariate_codes_preserve_probability_and_transport_sign():
    probability = np.array([[0.0, 1.0, 1.0]])
    transport = np.array([[0.0, -2.0, 2.0]])
    codes, normalized_probability, normalized_signed = signed_bivariate_codes(
        probability, transport
    )
    assert codes.shape == probability.shape
    assert normalized_probability[0, 1] == normalized_probability[0, 2]
    assert normalized_signed[0, 1] == -normalized_signed[0, 2]
    assert normalized_signed[0, 1] < 0.0 < normalized_signed[0, 2]


def test_signed_bivariate_colors_distinguish_transport_direction():
    positive = signed_bivariate_reference_rgb(0.4, 0.8)
    negative = signed_bivariate_reference_rgb(0.4, -0.8)
    assert positive[2] > positive[0]
    assert negative[0] > negative[2]
    assert positive != negative
    assert len(signed_bivariate_reference_colorscale()) == (
        BIVARIATE_LEVELS * (2 * BIVARIATE_LEVELS - 1)
    )
