import sympy as sp
from sympy import abc, Rational
from sympy.stats import Normal, density

import symbols as sym


# w equations
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

def w_bar(alpha=sp.abc.alpha, delta=sp.abc.delta, w_1=sym.w_1, w_2=sym.w_2):
    return ((1 - delta) * alpha * w_1
            + (1 - delta) * (1 - alpha) * w_2
            + delta * (alpha * w_1 + (1 - alpha) * w_2))


def w_prime_2_bar(alpha=sp.abc.alpha, delta=sp.abc.delta, w_1=sym.w_1, w_2=sym.w_2,
                  w_bar=sym.w_bar, sigma_w=sym.sigma_w, sigma_w_3=sym.sigma_w_3):
    return (((1 - delta) * alpha * ((w_1 - w_bar) ** 2 + sigma_w ** 2)) +
            ((1 - delta) * (1 - alpha) * ((w_2 - w_bar) ** 2 + sigma_w ** 2)) +
            (delta * sigma_w_3 ** 2))


def w_prime_3_bar(alpha=sp.abc.alpha, delta=sp.abc.delta, w_1=sym.w_1, w_2=sym.w_2,
                  w_bar=sym.w_bar, sigma_w=sym.sigma_w):
    return (((1 - delta) * alpha * ((w_1 - w_bar) ** 3 +
                                    3 * sigma_w ** 2 * (w_1 - w_bar))) +
            ((1 - delta) * (1 - alpha) * ((w_2 - w_bar) ** 3 +
                                          3 * sigma_w ** 2 * (w_2 - w_bar))))


def w_prime_4_bar(w_prime_2_bar=sym.w_prime_2_bar,
                  w_prime_3_bar=sym.w_prime_3_bar,
                  delta=sp.abc.delta,
                  sigma_tilde_w=sym.sigma_tilde_w,
                  sigma_w_3=sym.sigma_w_3):
    return (w_prime_2_bar ** 2 *
            ((1 - delta * (sigma_w_3 ** 2 / w_prime_2_bar)) ** 2 / (1 - delta)) *
            (3 * sigma_tilde_w ** 4 +
             6 * (1 - sigma_tilde_w ** 2) *
             sigma_tilde_w ** 2 +
             (1 - sigma_tilde_w ** 2) ** 2) +
            ((1 / (1 - sigma_tilde_w ** 2)) *
             (1 / (1 - delta * (sigma_w_3 ** 2 / w_prime_2_bar))) *
             (w_prime_3_bar ** 2 / w_prime_2_bar)) +
            (delta * 3 * sigma_w_3 ** 4))


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# theta_l equations
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

def theta_l_bar(alpha=sp.abc.alpha, delta=sp.abc.delta,
                theta_l_1=sym.theta_l_1, theta_l_2=sym.theta_l_2):
    return ((1 - delta) * alpha * theta_l_1
            + (1 - delta) * (1 - alpha) * theta_l_2
            + delta * (alpha * theta_l_1 + (1 - alpha) * theta_l_2))


def theta_l_prime_2_bar(alpha=sp.abc.alpha,
                        delta=sp.abc.delta,
                        theta_l_1=sym.theta_l_1,
                        theta_l_2=sym.theta_l_2,
                        theta_l_bar=sym.theta_l_bar,
                        sigma_theta_l_1=sym.sigma_theta_l_1,
                        sigma_theta_l_2=sym.sigma_theta_l_2,
                        sigma_theta_l_3_l=sym.sigma_theta_l_3):
    return (((1 - delta) * alpha * ((theta_l_1 - theta_l_bar) ** 2 + sigma_theta_l_1 ** 2)) +
            ((1 - delta) * (1 - alpha) *
             ((theta_l_2 - theta_l_bar) ** 2 + sigma_theta_l_2 ** 2)) +
            (delta * sigma_theta_l_3_l ** 2))


def theta_l_prime_3_bar(delta=sp.abc.delta,
                        alpha=sp.abc.alpha,
                        theta_l_1=sym.theta_l_1,
                        theta_l_2=sym.theta_l_2,
                        theta_l_bar=sym.theta_l_bar,
                        sigma_theta_l_1=sym.sigma_theta_l_1,
                        sigma_theta_l_2=sym.sigma_theta_l_2):
    return (((1 - delta) * alpha *
             ((theta_l_1 - theta_l_bar) ** 3 +
              3 * sigma_theta_l_1 ** 2 * (theta_l_1 - theta_l_bar))) +
            ((1 - delta) * (1 - alpha) *
             ((theta_l_2 - theta_l_bar) ** 3 +
              3 * sigma_theta_l_2 ** 2 * (theta_l_2 - theta_l_bar))))


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# r_t equations
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

def r_t_bar(alpha=sp.abc.alpha, delta=sp.abc.delta, r_t_1=sym.r_t_1, r_t_2=sym.r_t_2):
    return ((1 - delta) * alpha * r_t_1
            + (1 - delta) * (1 - alpha) * r_t_2
            + delta * (alpha * r_t_1 + (1 - alpha) * r_t_2))


def r_t_prime_2_bar(delta=sp.abc.delta,
                    alpha=sp.abc.alpha,
                    r_t_1=sym.r_t_1,
                    r_t_2=sym.r_t_2,
                    r_t_bar=sym.r_t_bar,
                    sigma_r_t_1=sym.sigma_r_t_1,
                    sigma_r_t_2=sym.sigma_r_t_2,
                    sigma_r_t_3_t=sym.sigma_r_t_3):
    return (((1 - delta) * alpha * ((r_t_1 - r_t_bar) ** 2 + sigma_r_t_1 ** 2)) +
            ((1 - delta) * (1 - alpha) * ((r_t_2 - r_t_bar) ** 2 + sigma_r_t_2 ** 2)) +
            (delta * sigma_r_t_3_t ** 2))


def r_t_prime_3_bar(alpha=sp.abc.alpha,
                    delta=sp.abc.delta,
                    r_t_1=sym.r_t_1,
                    r_t_2=sym.r_t_2,
                    r_t_bar=sym.r_t_bar,
                    sigma_r_t_1=sym.sigma_r_t_1,
                    sigma_r_t_2=sym.sigma_r_t_2):
    return (((1 - delta) * alpha *
             ((r_t_1 - r_t_bar) ** 3 +
              3 * sigma_r_t_1 ** 2 * (r_t_1 - r_t_bar))) +
            ((1 - delta) * (1 - alpha) *
             ((r_t_2 - r_t_bar) ** 3 +
              3 * sigma_r_t_2 ** 2 * (r_t_2 - r_t_bar))))


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# Mixed equations
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

def w_prime_theta_l_prime_bar(delta=sp.abc.delta,
                              alpha=sp.abc.alpha,
                              w_1=sym.w_1,
                              w_2=sym.w_2,
                              w_bar=sym.w_bar,
                              theta_l_1=sym.theta_l_1,
                              theta_l_2=sym.theta_l_2,
                              theta_l_bar=sym.theta_l_bar,
                              cov_lambda_w_theta=sym.rho_w_theta_l * sym.sigma_w_3 * sym.sigma_theta_l_3):
    return (((1 - delta) * alpha * ((w_1 - w_bar) * (theta_l_1 - theta_l_bar))) +
            ((1 - delta) * (1 - alpha) * ((w_2 - w_bar) * (theta_l_2 - theta_l_bar)))
            + delta * cov_lambda_w_theta)


def w_prime_r_t_prime_bar(alpha=sp.abc.alpha,
                          delta=sp.abc.delta,
                          w_1=sym.w_1,
                          w_2=sym.w_2,
                          w_bar=sym.w_bar,
                          r_t_1=sym.r_t_1,
                          r_t_2=sym.r_t_2,
                          r_t_bar=sym.r_t_bar,
                          cov_lambda_w_r=sym.rho_w_r_t * sym.sigma_w_3 * sym.sigma_r_t_3):
    return (((1 - delta) * alpha * ((w_1 - w_bar) * (r_t_1 - r_t_bar))) +
            ((1 - delta) * (1 - alpha) * ((w_2 - w_bar) * (r_t_2 - r_t_bar)))
            + delta * cov_lambda_w_r)


def w_prime_2_theta_l_prime_bar(sigma_tilde_w=sym.sigma_tilde_w,
                                delta=sp.abc.delta,
                                lambda_w_theta=sym.lambda_w_theta,
                                lambda_w=sym.lambda_w,
                                w_prime_3_bar=sym.w_prime_3_bar,
                                w_prime_2_bar=sym.w_prime_2_bar,
                                w_prime_theta_l_prime=sym.w_prime_theta_l_prime_bar):
    return ((1 / (1 - sigma_tilde_w ** 2)) *
            ((1 - delta * lambda_w_theta) / (1 - delta * lambda_w)) *
            (w_prime_3_bar / w_prime_2_bar) *
            w_prime_theta_l_prime)


def w_prime_theta_l_prime_2_bar(delta=sp.abc.delta,
                                lambda_w_theta=sym.lambda_w_theta,
                                lambda_w=sym.lambda_w,
                                sigma_tilde_w=sym.sigma_tilde_w,
                                w_prime_3_bar=sym.w_prime_3_bar,
                                w_prime_2_bar=sym.w_prime_2_bar,
                                w_prime_theta_l_prime_bar=sym.w_prime_theta_l_prime_bar,
                                theta_l_prime_3_bar=sym.theta_l_prime_3_bar):
    return (Rational(2, 3) *
            ((1 - delta * lambda_w_theta) ** 2 / (1 - delta * lambda_w) ** 2) *
            (1 / (1 - sigma_tilde_w ** 2) ** 2) *
            (w_prime_3_bar / w_prime_2_bar ** 2) *
            w_prime_theta_l_prime_bar ** 2 +
            Rational(1, 3) *
            ((1 - delta * lambda_w) / (1 - delta * lambda_w_theta)) *
            (1 - sigma_tilde_w ** 2) *
            ((w_prime_2_bar * theta_l_prime_3_bar) / w_prime_theta_l_prime_bar))


def w_prime_theta_l_prime_2_bar_beta(
        sigma_tilde_w=sym.sigma_tilde_w,
        delta=sp.abc.delta,
        lambda_theta=sym.lambda_theta,
        lambda_w=sym.lambda_w,
        w_prime_3_bar=sym.w_prime_3_bar,
        w_prime_2_bar=sym.w_prime_2_bar,
        beta=sp.abc.beta,
        theta_l_prime_2_bar=sym.theta_l_prime_2_bar,
        lambda_w_theta=sym.lambda_w_theta,
        w_prime_theta_l_prime_bar=sym.w_prime_theta_l_prime_bar):
    from sympy import Rational
    return (
            (1 / (1 - sigma_tilde_w ** 2)) *
            ((1 - delta * lambda_theta) / (1 - delta * lambda_w)) *
            (w_prime_3_bar / w_prime_2_bar) *
            (
                    Rational(1, 3) * beta * theta_l_prime_2_bar +
                    (
                            ((1 - Rational(1, 3) * beta) / (1 - sigma_tilde_w ** 2)) *
                            ((1 - delta * lambda_w_theta) ** 2 /
                             ((1 - delta * lambda_w) * (1 - delta * lambda_theta))) *
                            (w_prime_theta_l_prime_bar ** 2 / w_prime_2_bar)
                    )
            )
    )


def w_prime_r_t_prime_theta_l_prime_bar(delta=sp.abc.delta,
                                        alpha=sp.abc.alpha,
                                        w_1=sym.w_1,
                                        w_2=sym.w_2,
                                        w_bar=sym.w_bar,
                                        r_t_1=sym.r_t_1,
                                        r_t_2=sym.r_t_2,
                                        r_t_bar=sym.r_t_bar,
                                        theta_l_1=sym.theta_l_1,
                                        theta_l_2=sym.theta_l_2,
                                        theta_l_bar=sym.theta_l_bar,
                                        r_r_t_theta_l=sym.r_r_t_theta_l,
                                        sigma_r_t_1=sym.sigma_r_t_1,
                                        sigma_r_t_2=sym.sigma_r_t_2,
                                        sigma_theta_l_1=sym.sigma_theta_l_1,
                                        sigma_theta_l_2=sym.sigma_theta_l_2):
    return (((1 - delta) * alpha * (w_1 - w_bar) *
             (
                     (r_t_1 - r_t_bar) * (theta_l_1 - theta_l_bar) +
                     r_r_t_theta_l * sigma_r_t_1 * sigma_theta_l_1
             )) +
            ((1 - delta) * (1 - alpha) * (w_2 - w_bar) *
             (
                     (r_t_2 - r_t_bar) * (theta_l_2 - theta_l_bar) +
                     r_r_t_theta_l * sigma_r_t_2 * sigma_theta_l_2
             )))


def w_prime_r_t_prime_theta_l_prime_bar_beta(
        beta=sp.abc.beta,
        sigma_tilde_w=sym.sigma_tilde_w,
        delta=sp.abc.delta,
        lambda_theta_r=sym.lambda_theta_r,
        lambda_w=sym.lambda_w,
        r_t_prime_theta_l_prime_bar=sym.r_t_prime_theta_l_prime_bar,
        w_prime_3_bar=sym.w_prime_3_bar,
        w_prime_2_bar=sym.w_prime_2_bar,
        lambda_w_r=sym.lambda_w_r,
        lambda_w_theta=sym.lambda_w_theta,
        w_prime_r_t_prime_bar=sym.w_prime_r_t_prime_bar,
        w_prime_theta_l_prime_bar=sym.w_prime_theta_l_prime_bar):
    from sympy import Rational
    return (
            (
                    ((Rational(1, 3) * beta) / (1 - sigma_tilde_w ** 2)) *
                    ((1 - delta * lambda_theta_r) / (1 - delta * lambda_w)) *
                    r_t_prime_theta_l_prime_bar *
                    (w_prime_3_bar / w_prime_2_bar)
            ) +
            (
                    ((1 - Rational(1, 3) * beta) / (1 - sigma_tilde_w ** 2) ** 2) *
                    (((1 - delta * lambda_w_r) * (1 - delta * lambda_w_theta)) /
                     ((1 - delta * lambda_w) ** 2)) *
                    w_prime_r_t_prime_bar *
                    w_prime_theta_l_prime_bar *
                    (w_prime_3_bar / w_prime_2_bar ** 2)
            )
    )


def w_prime_r_t_prime_theta_l_prime_bar_E(
        E=sp.abc.E,
        sigma_tilde_w=sym.sigma_tilde_w,
        delta=sp.abc.delta,
        lambda_theta_r=sym.lambda_theta_r,
        lambda_w=sym.lambda_w,
        r_t_prime_theta_l_prime_bar=sym.r_t_prime_theta_l_prime_bar,
        w_prime_3_bar=sym.w_prime_3_bar,
        w_prime_2_bar=sym.w_prime_2_bar,
        lambda_w_r=sym.lambda_w_r,
        lambda_w_theta=sym.lambda_w_theta,
        w_prime_r_t_prime_bar=sym.w_prime_r_t_prime_bar,
        w_prime_theta_l_prime_bar=sym.w_prime_theta_l_prime_bar):
    from sympy import Rational
    return (
            (
                    ((Rational(1, 2) * E) / (1 - sigma_tilde_w ** 2)) *
                    ((1 - delta * lambda_theta_r) / (1 - delta * lambda_w)) *
                    r_t_prime_theta_l_prime_bar *
                    (w_prime_3_bar / w_prime_2_bar)
            ) +
            (
                    ((1 - Rational(1, 2) * E) / (1 - sigma_tilde_w ** 2) ** 2) *
                    (((1 - delta * lambda_w_r) * (1 - delta * lambda_w_theta)) /
                     ((1 - delta * lambda_w) ** 2)) *
                    w_prime_r_t_prime_bar *
                    w_prime_theta_l_prime_bar *
                    (w_prime_3_bar / w_prime_2_bar ** 2)
            )
    )


def r_t_prime_theta_l_prime_bar(
        alpha=sp.abc.alpha, delta=sp.abc.delta,
        r_t_1=sym.r_t_1, r_t_2=sym.r_t_2, r_t_prime_bar=sym.r_t_bar,
        theta_l_1=sym.theta_l_1, theta_l_2=sym.theta_l_2,
        theta_l_bar=sym.theta_l_bar,
        r_r_t_theta_l=sym.r_r_t_theta_l,
        sigma_r_t_1=sym.sigma_r_t_1, sigma_r_t_2=sym.sigma_r_t_2,
        sigma_theta_l_1=sym.sigma_theta_l_1,
        sigma_theta_l_2=sym.sigma_theta_l_2,
        cov_lambda_r_theta=sym.rho_theta_l_r_t * sym.sigma_theta_l_3 * sym.sigma_r_t_3):
    return ((1 - delta) * alpha * (
            (r_t_1 - r_t_prime_bar) * (theta_l_1 - theta_l_bar) +
            r_r_t_theta_l * sigma_r_t_1 * sigma_theta_l_1) +
            ((1 - delta) * (1 - alpha) * (
                    (r_t_2 - r_t_prime_bar) * (theta_l_2 - theta_l_bar) +
                    r_r_t_theta_l * sigma_r_t_2 * sigma_theta_l_2)) +
            delta * cov_lambda_r_theta)


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# Distributions
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

G_1_theta_l = Normal(name='G_1_theta_l', mean=sym.theta_l_1, std=sym.sigma_theta_l_1)
G_1_theta_l_density = density(G_1_theta_l)(sym.theta_l)

G_2_theta_l = Normal(name='G_2_theta_l', mean=sym.theta_l_2, std=sym.sigma_theta_l_2)
G_2_theta_l_density = density(G_2_theta_l)(sym.theta_l)

G_3_theta_l = Normal(name='G_3_theta_l', mean=sym.theta_l_bar, std=sym.sigma_theta_l_3)
G_3_theta_l_density = density(G_3_theta_l)(sym.theta_l)

G_theta = ((1 - sp.abc.delta) * sp.abc.alpha * G_1_theta_l_density +
           (1 - sp.abc.delta) * (1 - sp.abc.alpha) * G_2_theta_l_density +
           sp.abc.delta * G_3_theta_l_density)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

G_1_w = Normal(name='G_1_w', mean=sym.w_1, std=sym.sigma_w)
G_1_w_density = density(G_1_w)(sp.abc.w)

G_2_w = Normal(name='G_2_w', mean=sym.w_2, std=sym.sigma_w)
G_2_w_density = density(G_2_w)(sp.abc.w)

G_3_w = Normal(name='G_3_w', mean=sym.w_bar, std=sym.sigma_w_3)
G_3_w_density = density(G_3_w)(sp.abc.w)

G_w = ((1 - sp.abc.delta) * sp.abc.alpha * G_1_w_density +
       (1 - sp.abc.delta) * (1 - sp.abc.alpha) * G_2_w_density +
       sp.abc.delta * G_3_w_density)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

G_1_w_theta = Normal(name='G_1_w_theta', mean=sp.Matrix([sym.w_1, sym.theta_l_1]),
                     std=sp.Matrix([[sym.sigma_w ** 2, 0], [0, sym.sigma_theta_l_1 ** 2]]))
G_1_w_theta_density = density(G_1_w_theta)(sp.abc.w, sp.abc.theta)

G_2_w_theta = Normal(name='G_2_w_theta', mean=sp.Matrix([sym.w_2, sym.theta_l_2]),
                     std=sp.Matrix([[sym.sigma_w ** 2, 0], [0, sym.sigma_theta_l_2 ** 2]]))
G_2_w_theta_density = density(G_2_w_theta)(sp.abc.w, sp.abc.theta)

G_3_w_theta = Normal(name='G_3_w_theta', mean=sp.Matrix([sym.w_bar, sym.theta_l_bar]),
                     std=sp.Matrix([
                         [sym.sigma_w_3 ** 2,
                          sym.rho_w_theta_l * sym.sigma_w_3 * sym.sigma_theta_l_3],
                         [sym.rho_w_theta_l * sym.sigma_w_3 * sym.sigma_theta_l_3,
                          sym.sigma_theta_l_3 ** 2]
                     ]))
G_3_w_theta_density = sp.simplify(density(G_3_w_theta)(sp.abc.w, sp.abc.theta))

G_w_theta = ((1 - sp.abc.delta) * sp.abc.alpha * G_1_w_theta_density +
             (1 - sp.abc.delta) * (1 - sp.abc.alpha) * G_2_w_theta_density +
             sp.abc.delta * G_3_w_theta_density)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

mu_1_theta_l_r_t = sp.Matrix([sym.theta_l_1, sym.r_t_1])
Sigma_1_theta_l_r_t = sp.Matrix([[sym.sigma_theta_l_1 ** 2,
                                  sym.r_r_t_theta_l * sym.sigma_theta_l_1 * sym.sigma_r_t_1],
                                 [sym.r_r_t_theta_l * sym.sigma_theta_l_1 * sym.sigma_r_t_1,
                                  sym.sigma_r_t_1 ** 2]])

G_1_theta_l_r_t = Normal(name='G_1_theta_l_r_t',
                         mean=mu_1_theta_l_r_t,
                         std=Sigma_1_theta_l_r_t)

G_1_theta_l_r_t_density = density(G_1_theta_l_r_t)(sym.theta_l, sym.r_t)

mu_2_theta_l_r_t = sp.Matrix([sym.theta_l_2, sym.r_t_2])
Sigma_2_theta_l_r_t = sp.Matrix(
    [[sym.sigma_theta_l_2 ** 2,
      sym.r_r_t_theta_l * sym.sigma_theta_l_2 * sym.sigma_r_t_2],
     [sym.r_r_t_theta_l * sym.sigma_theta_l_2 * sym.sigma_r_t_2,
      sym.sigma_r_t_2 ** 2]])

G_2_theta_l_r_t = Normal(name='G_2_theta_l_r_t',
                         mean=mu_2_theta_l_r_t,
                         std=Sigma_2_theta_l_r_t)

G_2_theta_l_r_t_density = density(G_2_theta_l_r_t)(sym.theta_l, sym.r_t)

mu_3_theta_l_r_t = sp.Matrix([sym.theta_l_bar, sym.r_t_bar])
Sigma_3_theta_l_r_t = sp.Matrix(
    [[sym.sigma_theta_l_3 ** 2,
      sym.rho_theta_l_r_t * sym.sigma_theta_l_3 * sym.sigma_r_t_3],
     [
         sym.rho_theta_l_r_t * sym.sigma_theta_l_3 * sym.sigma_r_t_3,
         sym.sigma_r_t_3 ** 2]])

G_3_theta_l_r_t = Normal(name='G_3_theta_l_r_t',
                         mean=mu_3_theta_l_r_t,
                         std=Sigma_3_theta_l_r_t)

G_3_theta_l_r_t_density = density(G_3_theta_l_r_t)(sym.theta_l, sym.r_t)

G_theta_l_r_t = ((1 - sp.abc.delta) * sp.abc.alpha * G_1_theta_l_r_t_density +
                 (1 - sp.abc.delta) * (1 - sp.abc.alpha) * G_2_theta_l_r_t_density +
                 sp.abc.delta * G_3_theta_l_r_t_density)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

mu_1_w_theta_l_r_t = sp.Matrix([sym.w_1, sym.theta_l_1, sym.r_t_1])
Sigma_1_w_theta_l_r_t = sp.Matrix(
    [[sym.sigma_w ** 2,
      0,
      0],
     [0,
      sym.sigma_theta_l_1 ** 2,
      sym.r_r_t_theta_l * sym.sigma_theta_l_1 * sym.sigma_r_t_1],
     [0,
      sym.r_r_t_theta_l * sym.sigma_theta_l_1 * sym.sigma_r_t_1,
      sym.sigma_r_t_1 ** 2]])

G_1_w_theta_l_r_t = Normal(name='G_1_w_theta_l_r_t',
                           mean=mu_1_w_theta_l_r_t,
                           std=Sigma_1_w_theta_l_r_t)

G_1_w_theta_l_r_t_density = density(G_1_w_theta_l_r_t)(sp.abc.w, sym.theta_l, sym.r_t)

mu_2_w_theta_l_r_t = sp.Matrix([sym.w_2, sym.theta_l_2, sym.r_t_2])
Sigma_2_w_theta_l_r_t = sp.Matrix([[sym.sigma_w ** 2,
                                    0,
                                    0],
                                   [0,
                                    sym.sigma_theta_l_2 ** 2,
                                    sym.r_r_t_theta_l * sym.sigma_theta_l_2 * sym.sigma_r_t_2],
                                   [0,
                                    sym.r_r_t_theta_l * sym.sigma_theta_l_2 * sym.sigma_r_t_2,
                                    sym.sigma_r_t_2 ** 2]])

G_2_w_theta_l_r_t = Normal(name='G_2_w_theta_l_r_t',
                           mean=mu_2_w_theta_l_r_t,
                           std=Sigma_2_w_theta_l_r_t)

G_2_w_theta_l_r_t_density = density(G_2_w_theta_l_r_t)(sp.abc.w, sym.theta_l, sym.r_t)

mu_3_w_theta_l_r_t = sp.Matrix([sym.w_bar, sym.theta_l_bar, sym.r_t_bar])
Sigma_3_w_theta_l_r_t = sp.Matrix(
    [[sym.sigma_w_3 ** 2,
      sym.rho_w_theta_l * sym.sigma_w_3 * sym.sigma_theta_l_3,
      sym.rho_w_r_t * sym.sigma_theta_l_3 * sym.sigma_r_t_3],
     [
         sym.rho_w_theta_l * sym.sigma_w_3 * sym.sigma_theta_l_3,
         sym.sigma_theta_l_3 ** 2,
         sym.rho_theta_l_r_t * sym.sigma_w_3 * sym.sigma_r_t_3],
     [sym.rho_w_r_t * sym.sigma_theta_l_3 * sym.sigma_r_t_3,
      sym.rho_theta_l_r_t * sym.sigma_w_3 * sym.sigma_r_t_3,
      sym.sigma_r_t_3 ** 2]])

G_3_w_theta_l_r_t = Normal(name='G_3_w_theta_l_r_t',
                           mean=mu_3_w_theta_l_r_t,
                           std=Sigma_2_w_theta_l_r_t)
G_3_w_theta_l_r_t_density = density(G_3_w_theta_l_r_t)(sp.abc.w, sym.theta_l, sym.r_t)

G_w_theta_l_r_t = ((1 - sp.abc.delta) * sp.abc.alpha * G_1_w_theta_l_r_t_density +
                   (1 - sp.abc.delta) * (1 - sp.abc.alpha) * G_2_w_theta_l_r_t_density +
                   sp.abc.delta * G_3_w_theta_l_r_t_density)


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# sigma equations
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

def sigma_tilde_w(sigma_w=sym.sigma_w, w_prime_2_bar=sym.w_prime_2_bar,
                  delta=sp.abc.delta, lambda_w=sym.lambda_w):
    from sympy import sqrt
    return ((sigma_w / sqrt(w_prime_2_bar)) *
            (1 / sqrt((1 - delta * lambda_w) / (1 - delta))))


def sigma_tilde_r_t_1(sigma_r_t_1=sym.sigma_r_t_1, r_t_prime_2_bar=sym.r_t_prime_2_bar,
                      delta=sp.abc.delta, lambda_r=sym.lambda_r):
    from sympy import sqrt
    return ((sigma_r_t_1 / sqrt(r_t_prime_2_bar)) *
            (1 / sqrt((1 - delta * lambda_r) / (1 - delta))))


def sigma_tilde_r_t_2(sigma_r_t_2=sym.sigma_r_t_2, r_t_prime_2_bar=sym.r_t_prime_2_bar,
                      delta=sp.abc.delta, lambda_r=sym.lambda_r):
    from sympy import sqrt
    return ((sigma_r_t_2 / sqrt(r_t_prime_2_bar)) *
            (1 / sqrt((1 - delta * lambda_r) / (1 - delta))))


def sigma_tilde_theta_l_1(sigma_theta_l_1=sym.sigma_theta_l_1,
                          theta_l_prime_2_bar=sym.theta_l_prime_2_bar,
                          delta=sp.abc.delta, lambda_theta=sym.lambda_theta):
    from sympy import sqrt
    return ((sigma_theta_l_1 / sqrt(theta_l_prime_2_bar)) *
            (1 / sqrt((1 - delta * lambda_theta) / (1 - delta))))


def sigma_tilde_theta_l_2(sigma_theta_l_2=sym.sigma_theta_l_2,
                          theta_l_prime_2_bar=sym.theta_l_prime_2_bar,
                          delta=sp.abc.delta, lambda_theta=sym.lambda_theta):
    from sympy import sqrt
    return ((sigma_theta_l_2 / sqrt(theta_l_prime_2_bar)) *
            (1 / sqrt((1 - delta * lambda_theta) / (1 - delta))))


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# lambda equations
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

def lambda_w(sigma_w_3=sym.sigma_w_3, w_prime_2_bar=sym.w_prime_2_bar):
    return sigma_w_3 ** 2 / w_prime_2_bar


def lambda_theta(sigma_theta_l_3=sym.sigma_theta_l_3,
                 theta_l_prime_2_bar=sym.theta_l_prime_2_bar):
    return sigma_theta_l_3 ** 2 / theta_l_prime_2_bar


def lambda_r(sigma_r_t_3=sym.sigma_r_t_3,
             r_t_prime_2_bar=sym.r_t_prime_2_bar):
    return sigma_r_t_3 ** 2 / r_t_prime_2_bar


def lambda_w_theta(cov_lambda_w_theta=sym.rho_w_theta_l * sym.sigma_w_3 * sym.sigma_theta_l_3,
                   w_prime_theta_l_prime_bar=sym.w_prime_theta_l_prime_bar):
    return cov_lambda_w_theta / w_prime_theta_l_prime_bar


def lambda_w_r(cov_lambda_w_r=sym.rho_w_r_t * sym.sigma_w_3 * sym.sigma_r_t_3,
               w_prime_r_t_prime_bar=sym.w_prime_r_t_prime_bar):
    return cov_lambda_w_r / w_prime_r_t_prime_bar


def lambda_r_theta(cov_lambda_r_theta=sym.rho_theta_l_r_t * sym.sigma_r_t_3 * sym.sigma_theta_l_3,
                   r_t_prime_theta_l_prime_bar=sym.r_t_prime_theta_l_prime_bar):
    return cov_lambda_r_theta / r_t_prime_theta_l_prime_bar


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# sk
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

def sk_theta_l_hat_beta(sk_w_hat=sym.sk_w_hat, c_hat_w_theta_l=sym.c_w_theta_l_hat,
                        beta=sp.abc.beta):
    return sk_w_hat * c_hat_w_theta_l * (beta + (1 - beta) * c_hat_w_theta_l ** 2)


def sk_theta_l_hat(theta_l_prime_3_bar=sym.theta_l_prime_3_bar,
                   theta_l_prime_2_bar=sym.theta_l_prime_2_bar,
                   delta=sp.abc.delta,
                   lambda_theta=sym.lambda_theta):
    from sympy import Rational
    return ((theta_l_prime_3_bar / theta_l_prime_2_bar ** Rational(3, 2)) *
            (1 / ((1 - delta * lambda_theta) / (1 - delta)) ** Rational(3, 2)) *
            (1 / (1 - delta)))


def sk_r_t_hat(r_t_prime_3_bar=sym.r_t_prime_3_bar,
               r_t_prime_2_bar=sym.r_t_prime_2_bar,
               delta=sp.abc.delta,
               lambda_r=sym.lambda_r):
    from sympy import Rational
    return ((r_t_prime_3_bar / r_t_prime_2_bar ** Rational(3, 2)) *
            (1 / ((1 - delta * lambda_r) / (1 - delta)) ** Rational(3, 2)) *
            (1 / (1 - delta)))


def sk_w_hat(sigma_tilde_w=sym.sigma_tilde_w,
             w_prime_3_bar=sym.w_prime_3_bar,
             w_prime_2_bar=sym.w_prime_2_bar,
             delta=sp.abc.delta,
             lambda_w=sym.lambda_w):
    from sympy import Rational
    return ((1 / (1 - sigma_tilde_w ** 2) ** Rational(3, 2)) *
            (w_prime_3_bar / (w_prime_2_bar ** Rational(3, 2))) *
            (1 / ((1 - delta * lambda_w) / (1 - delta)) ** Rational(3, 2)) *
            (1 / (1 - delta)))


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# c

def c_w_theta_l_hat(sigma_w_tilde=sym.sigma_tilde_w,
                    w_prime_theta_l_prime_bar=sym.w_prime_theta_l_prime_bar,
                    w_prime_2_bar=sym.w_prime_2_bar,
                    theta_l_prime_2_bar=sym.theta_l_prime_2_bar,
                    delta=sp.abc.delta,
                    lambda_w=sym.lambda_w,
                    lambda_theta=sym.lambda_theta,
                    lambda_w_theta=sym.lambda_w_theta):
    from sympy import sqrt
    return ((1 / sqrt(1 - sigma_w_tilde ** 2)) *
            (w_prime_theta_l_prime_bar / ((sqrt(w_prime_2_bar)) * sqrt(theta_l_prime_2_bar))) *
            (1 / sqrt((1 - delta * lambda_w) / (1 - delta))) *
            (1 / sqrt((1 - delta * lambda_theta) / (1 - delta))) *
            ((1 - delta * lambda_w_theta) / (1 - delta)))


def c_w_r_t_hat(sigma_w_tilde=sym.sigma_tilde_w,
                w_prime_r_t_prime_bar=sym.w_prime_r_t_prime_bar,
                w_prime_2_bar=sym.w_prime_2_bar,
                r_t_prime_2_bar=sym.r_t_prime_2_bar,
                delta=sp.abc.delta,
                lambda_w=sym.lambda_w,
                lambda_r=sym.lambda_r,
                lambda_w_r=sym.lambda_w_r):
    from sympy import sqrt
    return ((1 / sqrt(1 - sigma_w_tilde ** 2)) *
            (w_prime_r_t_prime_bar / ((sqrt(w_prime_2_bar)) * sqrt(r_t_prime_2_bar))) *
            (1 / sqrt((1 - delta * lambda_w) / (1 - delta))) *
            (1 / sqrt((1 - delta * lambda_r) / (1 - delta))) *
            ((1 - delta * lambda_w_r) / (1 - delta)))


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# beta
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

def beta_w_theta_l(c_w_theta_l_hat=sym.c_w_theta_l_hat,
                   sk_w_hat=sym.sk_w_hat,
                   sk_theta_l_hat=sym.sk_theta_l_hat):
    return (((c_w_theta_l_hat ** 3 * sk_w_hat) - sk_theta_l_hat) /
            (c_w_theta_l_hat * sk_w_hat * (c_w_theta_l_hat ** 2 - 1)))


def beta_w_r_t(c_w_r_t_hat=sym.c_w_r_t_hat,
               sk_w_hat=sym.sk_w_hat,
               sk_r_t_hat=sym.sk_r_t_hat):
    return (((c_w_r_t_hat ** 3 * sk_w_hat) - sk_r_t_hat) /
            (c_w_r_t_hat * sk_w_hat * (c_w_r_t_hat ** 2 - 1)))


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# E
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

def E(alpha=sp.abc.alpha, xi=sp.abc.xi):
    from sympy import Rational
    return ((1 - Rational(1, 2) * ((2 * alpha) / (1 - 2 * alpha)) * xi) /
            (1 + Rational(1, 2) * xi))


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --


# xi
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

def xi(alpha=sp.abc.alpha, sigma_tilde_r_t_1=sym.sigma_tilde_r_t_1,
       sigma_tilde_r_t_2=sym.sigma_tilde_r_t_2, sigma_tilde_theta_l_1=sym.sigma_tilde_theta_l_1,
       sigma_tilde_theta_l_2=sym.sigma_tilde_theta_l_2):
    return (((1 - alpha) / alpha) *
            (sigma_tilde_r_t_1 / sigma_tilde_r_t_2) *
            (sigma_tilde_theta_l_1 / sigma_tilde_theta_l_2) - 1)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
