def get_permutation_df(w_prime_2_bar=True, w_prime_3_bar=True, sigma_tilde_w=True,
                       sk_w_hat=True, r_t_prime_2_bar=True, r_t_prime_3_bar=True, sk_r_t_hat=True,
                       theta_l_prime_2_bar=True, theta_l_prime_3_bar=True, sk_theta_l_hat=True,
                       w_prime_theta_l_prime_bar=True, w_prime_r_t_prime_bar=True,
                       r_t_prime_theta_l_prime_bar=True, c_w_theta_l_hat=True, c_w_r_t_hat=True,
                       beta_w_theta_l=True, beta_w_r_t=True):
    import symbols as sym
    import pandas as pd
    from itertools import product
    import sympy as sp
    from sympy import abc
    import checked_functions as c_f

    df = pd.DataFrame(
        product(
            [.4, .6],  # alpha
            [0, .5, .9],  # delta
            [3, 5],  # w_1
            [-1, -2],  # w_2
            [2, 6],  # r_t_1
            [-2, 1],  # r_t_2
            [4, 7],  # theta_l_1
            [0, 2],  # theta_l_2
            [.5],  # r_r_t_theta_l
            [1.35],  # sigma_w
            [1, 1.2],  # sigma_r_t_1
            [1, 1.2],  # sigma_r_t_2
            [1, 1.5],  # sigma_theta_l_1
            [1, 1.5],  # sigma_theta_l_2
            [.6],  # rho_w_r_t
            [.7],  # rho_w_theta_l
            [.8],  # rho_theta_l_r_t
            [.6],  # sigma_lambda_w
            [.7],  # sigma_lambda_theta_l
            [.8]  # sigma_lambda_r_t
        ),
        columns=[sp.abc.alpha,
                 sp.abc.delta,
                 sym.w_1,
                 sym.w_2,
                 sym.r_t_1,
                 sym.r_t_2,
                 sym.theta_l_1,
                 sym.theta_l_2,
                 sym.r_r_t_theta_l,
                 sym.sigma_w,
                 sym.sigma_r_t_1,
                 sym.sigma_r_t_2,
                 sym.sigma_theta_l_1,
                 sym.sigma_theta_l_2,
                 sym.rho_w_r_t,
                 sym.rho_w_theta_l,
                 sym.rho_theta_l_r_t,
                 sym.sigma_w_3,
                 sym.sigma_theta_l_3,
                 sym.sigma_r_t_3])

    if w_prime_2_bar:
        w_prime_2_bar = sp.lambdify([
            sp.abc.alpha,
            sp.abc.delta,
            sym.w_1,
            sym.w_2,
            sym.sigma_w,
            sym.sigma_w_3],
            c_f.w_prime_2_bar().subs({
                sym.w_bar: c_f.w_bar()
            }))

        df[sym.w_prime_2_bar] = (
            w_prime_2_bar(df[sp.abc.alpha],
                          df[sp.abc.delta],
                          df[sym.w_1],
                          df[sym.w_2],
                          df[sym.sigma_w],
                          df[sym.sigma_w_3]))

    if w_prime_3_bar:
        w_prime_3_bar = sp.lambdify([
            sp.abc.alpha,
            sp.abc.delta,
            sym.w_1,
            sym.w_2,
            sym.sigma_w],
            c_f.w_prime_3_bar().subs({
                sym.w_bar: c_f.w_bar()
            }))

        df[sym.w_prime_3_bar] = (
            w_prime_3_bar(df[sp.abc.alpha],
                          df[sp.abc.delta],
                          df[sym.w_1],
                          df[sym.w_2],
                          df[sym.sigma_w]))

    if sigma_tilde_w:
        if not w_prime_2_bar:
            raise KeyError('sigma_tilde_w needs w_prime_2_bar!')

        sigma_tilde_w = sp.lambdify([
            sym.w_prime_2_bar,
            sp.abc.delta,
            sym.sigma_w,
            sym.sigma_w_3],
            c_f.sigma_tilde_w().subs({
                sym.lambda_w: c_f.lambda_w()
            }))

        df[sym.sigma_tilde_w] = (
            sigma_tilde_w(df[sym.w_prime_2_bar],
                          df[sp.abc.delta],
                          df[sym.sigma_w],
                          df[sym.sigma_w_3]))

    if sk_w_hat:
        if not (sigma_tilde_w and w_prime_3_bar and w_prime_2_bar):
            raise KeyError('sk_w_hat needs sigma_tilde_w, w_prime_3_bar and w_prime_2_bar!')

        sk_w_hat = sp.lambdify([
            sym.sigma_tilde_w,
            sym.w_prime_3_bar,
            sym.w_prime_2_bar,
            sp.abc.delta,
            sym.sigma_w_3],
            c_f.sk_w_hat().subs({
                sym.lambda_w: c_f.lambda_w(),
            }))

        df[sym.sk_w_hat] = sk_w_hat(
            df[sym.sigma_tilde_w],
            df[sym.w_prime_3_bar],
            df[sym.w_prime_2_bar],
            df[sym.sigma_w_3],
            df[sym.sigma_tilde_w])

    if r_t_prime_2_bar:
        r_t_prime_2_bar = sp.lambdify([
            sp.abc.alpha,
            sp.abc.delta,
            sym.r_t_1,
            sym.r_t_2,
            sym.sigma_r_t_1,
            sym.sigma_r_t_2,
            sym.sigma_r_t_3],
            c_f.r_t_prime_2_bar().subs({
                sym.r_t_bar: c_f.r_t_bar()
            }))

        df[sym.r_t_prime_2_bar] = r_t_prime_2_bar(
            df[sp.abc.alpha],
            df[sp.abc.delta],
            df[sym.r_t_1],
            df[sym.r_t_2],
            df[sym.sigma_r_t_1],
            df[sym.sigma_r_t_2],
            df[sym.sigma_r_t_3])

    if r_t_prime_3_bar:
        r_t_prime_3_bar = sp.lambdify([
            sp.abc.alpha,
            sp.abc.delta,
            sym.r_t_1,
            sym.r_t_2,
            sym.sigma_r_t_1,
            sym.sigma_r_t_2],
            c_f.r_t_prime_3_bar().subs({
                sym.r_t_bar: c_f.r_t_bar(),
            }))

        df[sym.r_t_prime_3_bar] = r_t_prime_3_bar(
            df[sp.abc.alpha],
            df[sp.abc.delta],
            df[sym.r_t_1],
            df[sym.r_t_2],
            df[sym.sigma_r_t_1],
            df[sym.sigma_r_t_2])

    if sk_r_t_hat:
        if not (r_t_prime_3_bar and r_t_prime_2_bar):
            raise KeyError('sk_r_t_hat needs r_t_prime_3_bar r_t_prime_2_bar!')

        sk_r_t_hat = sp.lambdify([
            sym.r_t_prime_3_bar,
            sym.r_t_prime_2_bar,
            sp.abc.delta,
            sym.sigma_r_t_3],
            c_f.sk_r_t_hat().subs({
                sym.lambda_r: c_f.lambda_r(),
            }))

        df[sym.sk_r_t_hat] = sk_r_t_hat(
            df[sym.r_t_prime_3_bar],
            df[sym.r_t_prime_2_bar],
            df[sp.abc.delta],
            df[sym.sigma_r_t_3])

    if theta_l_prime_2_bar:
        theta_l_prime_2_bar = sp.lambdify([
            sp.abc.alpha,
            sp.abc.delta,
            sym.theta_l_1,
            sym.theta_l_2,
            sym.sigma_theta_l_1,
            sym.sigma_theta_l_2,
            sym.sigma_theta_l_3],
            c_f.theta_l_prime_2_bar().subs({
                sym.theta_l_bar: c_f.theta_l_bar()
            }))

        df[sym.theta_l_prime_2_bar] = theta_l_prime_2_bar(
            df[sp.abc.alpha],
            df[sp.abc.delta],
            df[sym.theta_l_1],
            df[sym.theta_l_2],
            df[sym.sigma_theta_l_1],
            df[sym.sigma_theta_l_2],
            df[sym.sigma_theta_l_3])

    if theta_l_prime_3_bar:
        theta_l_prime_3_bar = sp.lambdify([
            sp.abc.alpha,
            sp.abc.delta,
            sym.theta_l_1,
            sym.theta_l_2,
            sym.sigma_theta_l_1,
            sym.sigma_theta_l_2],
            c_f.theta_l_prime_3_bar().subs({
                sym.theta_l_bar: c_f.theta_l_bar(),
            }))

        df[sym.theta_l_prime_3_bar] = theta_l_prime_3_bar(
            df[sp.abc.alpha],
            df[sp.abc.delta],
            df[sym.theta_l_1],
            df[sym.theta_l_2],
            df[sym.sigma_theta_l_1],
            df[sym.sigma_theta_l_2])

    if sk_theta_l_hat:
        if not (theta_l_prime_3_bar and theta_l_prime_2_bar):
            raise KeyError('sk_theta_l_hat needs theta_l_prime_3_bar theta_l_prime_2_bar!')

        sk_theta_l_hat = sp.lambdify([
            sym.theta_l_prime_3_bar,
            sym.theta_l_prime_2_bar,
            sp.abc.delta,
            sym.sigma_theta_l_3],
            c_f.sk_theta_l_hat().subs({
                sym.lambda_theta: c_f.lambda_theta(),
            }))

        df[sym.sk_theta_l_hat] = sk_theta_l_hat(
            df[sym.theta_l_prime_3_bar],
            df[sym.theta_l_prime_2_bar],
            df[sp.abc.delta],
            df[sym.sigma_theta_l_3])

    if w_prime_theta_l_prime_bar:
        w_prime_theta_l_prime_bar = sp.lambdify([
            sp.abc.alpha,
            sp.abc.delta,
            sym.w_1,
            sym.w_2,
            sym.theta_l_1,
            sym.theta_l_2,
            sym.rho_w_theta_l,
            sym.sigma_w_3,
            sym.sigma_theta_l_3],
            c_f.w_prime_theta_l_prime_bar().subs({
                sym.w_bar: c_f.w_bar(),
                sym.theta_l_bar: c_f.theta_l_bar()
            }))

        df[sym.w_prime_theta_l_prime_bar] = w_prime_theta_l_prime_bar(
            df[sp.abc.alpha],
            df[sp.abc.delta],
            df[sym.w_1],
            df[sym.w_2],
            df[sym.theta_l_1],
            df[sym.theta_l_2],
            df[sym.rho_w_theta_l],
            df[sym.sigma_w_3],
            df[sym.sigma_theta_l_3])

    if w_prime_r_t_prime_bar:
        w_prime_r_t_prime_bar = sp.lambdify([
            sp.abc.alpha,
            sp.abc.delta,
            sym.w_1,
            sym.w_2,
            sym.r_t_1,
            sym.r_t_2,
            sym.rho_w_r_t,
            sym.sigma_w_3,
            sym.sigma_r_t_3],
            c_f.w_prime_r_t_prime_bar().subs({
                sym.w_bar: c_f.w_bar(),
                sym.r_t_bar: c_f.r_t_bar()
            }))

        df[sym.w_prime_r_t_prime_bar] = w_prime_r_t_prime_bar(
            df[sp.abc.alpha],
            df[sp.abc.delta],
            df[sym.w_1],
            df[sym.w_2],
            df[sym.r_t_1],
            df[sym.r_t_2],
            df[sym.rho_w_r_t],
            df[sym.sigma_w_3],
            df[sym.sigma_r_t_3])

    if r_t_prime_theta_l_prime_bar:
        r_t_prime_theta_l_prime_bar = sp.lambdify([
            sp.abc.alpha,
            sp.abc.delta,
            sym.r_t_1,
            sym.r_t_2,
            sym.theta_l_1,
            sym.theta_l_2,
            sym.r_r_t_theta_l,
            sym.sigma_r_t_1,
            sym.sigma_r_t_2,
            sym.sigma_theta_l_1,
            sym.sigma_theta_l_2,
            sym.rho_theta_l_r_t,
            sym.sigma_theta_l_3,
            sym.sigma_r_t_3
        ], c_f.r_t_prime_theta_l_prime_bar().subs({
            sym.r_t_bar: c_f.r_t_bar(),
            sym.theta_l_bar: c_f.theta_l_bar()
        }))

        df[sym.r_t_prime_theta_l_prime_bar] = r_t_prime_theta_l_prime_bar(
            df[sp.abc.alpha],
            df[sp.abc.delta],
            df[sym.r_t_1],
            df[sym.r_t_2],
            df[sym.theta_l_1],
            df[sym.theta_l_2],
            df[sym.r_r_t_theta_l],
            df[sym.sigma_r_t_1],
            df[sym.sigma_r_t_2],
            df[sym.sigma_theta_l_1],
            df[sym.sigma_theta_l_2],
            df[sym.rho_theta_l_r_t],
            df[sym.sigma_theta_l_3],
            df[sym.sigma_r_t_3])

    if c_w_theta_l_hat:
        if not (sigma_tilde_w and w_prime_r_t_prime_bar and w_prime_2_bar and theta_l_prime_2_bar):
            raise KeyError(
                'c_w_theta_l_hat needs sigma_tilde_w, w_prime_r_t_prime_bar, w_prime_2_bar '
                'and theta_l_prime_2_bar!')

        c_w_theta_l_hat = sp.lambdify([
            sym.sigma_tilde_w,
            sym.w_prime_theta_l_prime_bar,
            sym.w_prime_2_bar,
            sym.theta_l_prime_2_bar,
            sp.abc.delta,
            sym.sigma_w_3,
            sym.sigma_theta_l_3,
            sym.rho_w_theta_l],
            c_f.c_w_theta_l_hat().subs({
                sym.lambda_w: c_f.lambda_w(),
                sym.lambda_theta: c_f.lambda_theta(),
                sym.lambda_w_theta: c_f.lambda_w_theta()
            }))

        df[sym.c_w_theta_l_hat] = c_w_theta_l_hat(
            df[sym.sigma_tilde_w],
            df[sym.w_prime_theta_l_prime_bar],
            df[sym.w_prime_2_bar],
            df[sym.theta_l_prime_2_bar],
            df[sp.abc.delta],
            df[sym.sigma_w_3],
            df[sym.sigma_theta_l_3],
            df[sym.rho_w_theta_l])

    if c_w_r_t_hat:
        if not (sigma_tilde_w and w_prime_r_t_prime_bar and w_prime_2_bar and r_t_prime_2_bar):
            raise KeyError(
                'c_w_r_t_hat needs sigma_tilde_w, w_prime_r_t_prime_bar, w_prime_2_bar '
                'and r_t_prime_2_bar!')

        c_w_r_t_hat = sp.lambdify([
            sym.sigma_tilde_w,
            sym.w_prime_r_t_prime_bar,
            sym.w_prime_2_bar,
            sym.r_t_prime_2_bar,
            sp.abc.delta,
            sym.sigma_w_3,
            sym.sigma_r_t_3,
            sym.rho_w_r_t],
            c_f.c_w_r_t_hat().subs({
                sym.lambda_w: c_f.lambda_w(),
                sym.lambda_r: c_f.lambda_r(),
                sym.lambda_w_r: c_f.lambda_w_r()
            }))

        df[sym.c_w_r_t_hat] = c_w_r_t_hat(
            df[sym.sigma_tilde_w],
            df[sym.w_prime_r_t_prime_bar],
            df[sym.w_prime_2_bar],
            df[sym.r_t_prime_2_bar],
            df[sp.abc.delta],
            df[sym.sigma_w_3],
            df[sym.sigma_r_t_3],
            df[sym.rho_w_r_t])

    if beta_w_theta_l:
        if not (sk_theta_l_hat and sk_w_hat and c_w_theta_l_hat):
            raise KeyError(
                'beta_w_theta_l needs sk_theta_l_hat, sk_w_hat and c_w_theta_l_hat!')
        beta_w_theta_l = sp.lambdify([
            sym.sk_theta_l_hat,
            sym.sk_w_hat,
            sym.c_w_theta_l_hat],
            c_f.beta_w_theta_l())

        df[sym.beta_theta_l] = beta_w_theta_l(
            df[sym.sk_theta_l_hat],
            df[sym.sk_w_hat],
            df[sym.c_w_theta_l_hat])

    if beta_w_r_t:
        if not (sk_theta_l_hat and sk_r_t_hat and c_w_r_t_hat):
            raise KeyError(
                'beta_w_theta_l needs sk_r_t_hat, sk_w_hat and c_w_r_t_hat!')
        beta_w_r_t = sp.lambdify([
            sym.sk_r_t_hat,
            sym.sk_w_hat,
            sym.c_w_r_t_hat],
            c_f.beta_w_r_t())

        df[sym.beta_r_t] = beta_w_r_t(
            df[sym.sk_r_t_hat],
            df[sym.sk_w_hat],
            df[sym.c_w_r_t_hat])

    return df


def print_diagnostics(df):
    import numpy as np
    import sympy as sp
    from sympy import abc
    import symbols as sym

    def print_mean_error(diff_col_name='diff', delta_zero=False, sigma_r_t_same=False,
                         sigma_theta_l_same=False):

        if delta_zero:
            if sigma_r_t_same and sigma_theta_l_same:
                print('-', sp.abc.delta, '=', 0, 'and', sym.sigma_theta_l_1, '=',
                      sym.sigma_theta_l_2, '=', sym.sigma_r_t_1, '=', sym.sigma_r_t_2, ':',
                      np.mean(df.loc[(df[sp.abc.delta] == 0) &
                                     (df[sym.sigma_r_t_1] == df[sym.sigma_r_t_2]) &
                                     (df[sym.sigma_theta_l_1] == df[
                                         sym.sigma_theta_l_2]), diff_col_name]),
                      '\n')
                return

            if sigma_r_t_same:
                print('-', sp.abc.delta, '=', 0, 'and', sym.sigma_r_t_1, '=', sym.sigma_r_t_2, ':',
                      np.mean(df.loc[(df[sp.abc.delta] == 0) &
                                     (df[sym.sigma_r_t_1] == df[sym.sigma_r_t_2]), diff_col_name]),
                      '\n')
                return

            if sigma_theta_l_same:
                print('-', sp.abc.delta, '=', 0, 'and', sym.sigma_theta_l_1, '=',
                      sym.sigma_theta_l_2, ':',
                      np.mean(df.loc[(df[sp.abc.delta] == 0) &
                                     (df[sym.sigma_theta_l_1] == df[sym.sigma_theta_l_2]), diff_col_name]),
                      '\n')
                return

            print('-', sp.abc.delta, '=', 0, ':',
                  np.mean(df.loc[(df[sp.abc.delta] == 0), diff_col_name]),
                  '\n')
            return

        else:
            if sigma_r_t_same and sigma_theta_l_same:
                print('-', sym.sigma_theta_l_1, '=', sym.sigma_theta_l_2, '=',
                      sym.sigma_r_t_1, '=', sym.sigma_r_t_2, ':',
                      np.mean(df.loc[(df[sym.sigma_r_t_1] == df[sym.sigma_r_t_2]) &
                                     (df[sym.sigma_theta_l_1] == df[
                                         sym.sigma_theta_l_2]), diff_col_name]),
                      '\n')
                return

            if sigma_r_t_same:
                print('-', sym.sigma_r_t_1, '=', sym.sigma_r_t_2, ':',
                      np.mean(df.loc[(df[sym.sigma_r_t_1] == df[sym.sigma_r_t_2]), diff_col_name]),
                      '\n')
                return

            if sigma_theta_l_same:
                print('-', sym.sigma_theta_l_1, '=', sym.sigma_theta_l_2, ':',
                      np.mean(df.loc[(df[sym.sigma_theta_l_1] == df[sym.sigma_theta_l_2]), diff_col_name]),
                      '\n')
                return

            print('- overall:', np.mean(df[diff_col_name]), '\n')

    def print_percentage_error(diff_col_name='diff', greater_as=1e-14):
        print('Percentage of values greater than', greater_as, ':\n',
              (len(df[abs(df[diff_col_name]) > greater_as]) / len(df)) * 100, '%\n')

    def print_all_error_perms(diff_col_name='diff'):
        print_mean_error(diff_col_name=diff_col_name)
        print_mean_error(diff_col_name=diff_col_name, delta_zero=True)
        print_mean_error(diff_col_name=diff_col_name, delta_zero=True, sigma_r_t_same=True)
        print_mean_error(diff_col_name=diff_col_name, delta_zero=True, sigma_theta_l_same=True)
        print_mean_error(diff_col_name=diff_col_name, delta_zero=True,
                         sigma_r_t_same=True, sigma_theta_l_same=True)
        print_mean_error(diff_col_name=diff_col_name, sigma_r_t_same=True)
        print_mean_error(diff_col_name=diff_col_name, sigma_theta_l_same=True)
        print_mean_error(diff_col_name=diff_col_name, sigma_r_t_same=True, sigma_theta_l_same=True)

    if 'diff_theta' in df.columns:
        print('Mean errors between lhs and rhs in', sp.abc.theta, ':\n')
        print_all_error_perms(diff_col_name='diff_theta')
        print_percentage_error(diff_col_name='diff_theta')

    if 'diff_r' in df.columns:
        print('Mean errors between lhs and rhs in', sp.abc.r, ':\n')
        print_all_error_perms(diff_col_name='diff_r')
        print_percentage_error(diff_col_name='diff_r')

    if 'diff' in df.columns:
        print('Mean errors between lhs and rhs:\n')
        print_all_error_perms(diff_col_name='diff')
        print_percentage_error(diff_col_name='diff')
