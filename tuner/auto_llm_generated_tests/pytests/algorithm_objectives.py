"""Synthetic objectives from src/G_unit_test_types/tuner_tests.F90.

Goldstein–Price, Rastrigin, Himmelblau, Eggholder and Schaffer retain the
Fortran formulas and search domains. Rastrigin uses five coordinates, as in
that suite. The functions accept a physical parameter vector and return a
scalar objective; gradient=True also returns an independently written analytic
gradient. Value-only evaluations do not calculate derivatives. They need
no CLUBB case or JAX environment. Known minima define score gaps. Hit predicates
match the actual Fortran driver assertions;
Himmelblau accepts any point with sufficiently small loss, including its
distinct global minimizers. Eggholder's square-root cusps have no finite
derivative, so gradient evaluation rejects those exact points while scalar
evaluation remains available.

BENCHMARKS supplies defaults, bounds and known minimizers to the generated
Python algorithm tests. Gradients are verified separately using centered
finite differences at nonsingular interior points. These tests exercise the
Python fill/tell strategies. They use Fortran success rules, but do not claim
Fortran ESA parity or comparable 10,000-trial calibration.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np


def goldstein_price(point, *, gradient=False):
    """Return Goldstein–Price's loss and, when requested, its analytic partials."""
    x, y = point
    u = x + y + 1.
    v = 19. - 14.*x + 3.*x*x - 14.*y + 6.*x*y + 3.*y*y
    w = 2.*x - 3.*y
    z = 18. - 32.*x + 12.*x*x + 48.*y - 36.*x*y + 27.*y*y
    a, b = 1. + u*u*v, 30. + w*w*z
    if not gradient:
        return float(a*b)
    da = 2.*u*v*np.ones(2) + u*u*np.array([-14.+6.*x+6.*y, -14.+6.*x+6.*y])
    db = 2.*w*z*np.array([2., -3.]) + w*w*np.array([-32.+24.*x-36.*y, 48.-36.*x+54.*y])
    return float(a*b), da*b + a*db


def rastrigin(point, *, gradient=False):
    """Return Rastrigin's loss and optional oscillatory partials in n dimensions."""
    x = np.asarray(point)
    loss = float(10.*len(x) + np.sum(x*x - 10.*np.cos(2.*np.pi*x)))
    return (loss, 2.*x + 20.*np.pi*np.sin(2.*np.pi*x)) if gradient else loss


def himmelblau(point, *, gradient=False):
    """Return Himmelblau's loss and optional partials; four global minimizers exist."""
    x, y = point
    a, b = x*x + y - 11., x + y*y - 7.
    loss = float(a*a + b*b)
    return (loss, np.array([4.*x*a + 2.*b, 2.*a + 4.*y*b])) if gradient else loss


def eggholder(point, *, gradient=False):
    """Return Eggholder's loss; reject undefined cusp derivatives only if requested."""
    x, y = point
    u, v = x/2. + y + 47., x - y - 47.
    ru, rv = np.sqrt(abs(u)), np.sqrt(abs(v))
    su, sv = np.sin(ru), np.sin(rv)
    loss = float(-(y+47.)*su - x*sv)
    if not gradient:
        return loss
    if u == 0. or v == 0.:
        raise ValueError('Eggholder has no finite gradient at this square-root cusp')
    du, dv = np.cos(ru)*np.sign(u)/(2.*ru), np.cos(rv)*np.sign(v)/(2.*rv)
    return loss, np.array([-(y+47.)*du/2. - sv - x*dv,
                           -su - (y+47.)*du + x*dv])


def schaffer(point, *, gradient=False):
    """Return Schaffer's second loss and optional partials of its radial denominator."""
    x, y = point
    q = x*x - y*y
    a, b = np.sin(q)**2 - .5, 1. + .001*(x*x+y*y)
    loss = float(.5 + a/(b*b))
    if not gradient:
        return loss
    da = np.sin(2.*q)*np.array([2.*x, -2.*y])
    db = .002*np.array([x, y])
    return loss, da/(b*b) - 2.*a*db/(b*b*b)


def goldstein_price_hit(point, loss):
    """Recognize Goldstein–Price's minimum within inclusive coordinate tolerances."""
    return abs(point[0]) <= 1.e-2 and abs(point[1]+1.) <= 1.e-2


def rastrigin_hit(point, loss):
    """Recognize Rastrigin's origin only when all coordinates meet the strict tolerance."""
    return np.all(np.abs(point) < .1)


def himmelblau_hit(point, loss):
    """Accept any Himmelblau minimizer meeting Fortran's strict loss threshold."""
    return loss < 1.e-2


def eggholder_hit(point, loss):
    """Recognize Eggholder's boundary minimum using asymmetric strict tolerances."""
    return abs(point[0]-512.) < 1.e-2 and abs(point[1]-404.23) < 1.


def schaffer_hit(point, loss):
    """Recognize Schaffer's global target at or below the inclusive loss threshold."""
    return loss <= 1.e-5


@dataclass(frozen=True)
class Objective:
    """Describe a scalar/gradient callback, baseline and diagnostic global-hit rule."""

    name: str
    evaluate: Callable
    lower: float
    upper: float
    initial: tuple
    minimizer: tuple
    minimum: float
    hit_test: Callable
    hit_description: str


BENCHMARKS = (
    Objective('goldstein_price', goldstein_price, -2., 2., (1.4, 1.4), (0., -1.), 3.,
              goldstein_price_hit, '|x0| <= 0.01 and |x1 + 1| <= 0.01'),
    Objective('rastrigin', rastrigin, -5.12, 5.12, (3.7,)*5, (0.,)*5, 0.,
              rastrigin_hit, 'all |xi| < 0.1 (five coordinates)'),
    Objective('himmelblau', himmelblau, -5., 5., (0., 0.), (3., 2.), 0.,
              himmelblau_hit, 'f < 0.01 (any minimizer)'),
    Objective('eggholder', eggholder, -512., 512., (0., 0.), (512., 404.231805), -959.64066272,
              eggholder_hit, '|x0 - 512| < 0.01 and |x1 - 404.23| < 1'),
    Objective('schaffer', schaffer, -100., 100., (10., 10.), (0., 0.), 0.,
              schaffer_hit, 'f <= 1e-5'),
)
