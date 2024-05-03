# Adding a third Gaussian to CLUBB

In this folder there are several files which check the equations written in my masters thesis,
as well as in the document on overleaf.

## Helpers

- [checked_functions.py](checked_functions.py) is a document where all functions are defined.
  One usually imports this document and just uses the functions in there.
- [symbols.py](symbols.py) is the document where all the used symbols for computing things with ``sympy`` are defined. (These should also be imported into the calculations, e.g. jupyter notebooks.)

## Integral checks

All jupyter notebooks are basically checking if an integral on the "lhs" of an equation in numerically
or analytically equal to the "rhs" of an equation.

### Symbolic checks

We are checking whether the integrals equal the expressions in terms of pdf parameters.
The exception is $\overline{w'^4}$ which checks the rhs based on lower-order moments.
The symbolic checks are the following:

- [$\overline{\theta'_l^2}$](theta_l_prime_2_bar.ipynb)
- [$\overline{\theta'_l^3}$](theta_l_prime_3_bar.ipynb)
- [$\overline{w}$](w_bar.ipynb)
- [$\overline{w'^2}$](w_prime_2_bar.ipynb)
- [$\overline{w'^3}$](w_prime_3_bar.ipynb)
- [$\overline{w'^4}$](w_prime_4_bar.ipynb)

### Numeric checks

The numeric checks are the following:

- [$\overline{\theta_l'^3}\beta$](theta_l_prime_3_bar_beta.ipynb) (This one worked)
- [$\overline{w'r_t'\theta_l'}$](w_prime_r_t_prime_theta_l_prime_bar.ipynb) (This one is in terms of pdf parameters)
- [$\overline{w'r_t'\theta_l'}\beta$](w_prime_r_t_prime_theta_l_prime_bar_beta.ipynb) (That one did not work because of different $\beta$s)
- [$\overline{w'r_t'\theta_l'}E$](w_prime_r_t_prime_theta_l_prime_bar_E.ipynb) (That one did not work also)
- [$\overline{w'\theta_l'^2}$](w_prime_theta_l_prime_2_bar.ipynb) (This one is in terms of moments)
- [$\overline{w'^2\theta_l'}$](w_prime_2_theta_l_prime_bar.ipynb) (This one is in terms of moments)
- [$\overline{w'\theta_l'}$](w_prime_theta_l_prime_bar.ipynb) (This one is in terms of pdf parameters)

## How to use the documents

As said before, [checked_functions.py](checked_functions.py) is a document where all functions are defined.
The document [symbols.py](symbols.py) is the document where all the used symbols for computing things
with [``sympy``](https://www.sympy.org/en/index.html) are defined.
One usually calls a function with an empty list of arguments to get the function with the symbols already put in.
If you want to have a nice printout of the obtained equation, you should
import [``display``](https://ipython.readthedocs.io/en/stable/api/generated/IPython.display.html)
from ``IPython.display`` first.
Afterwards, one can use the function ``display(..)`` with a "rhs" and a "lhs".
At first, the "rhs" is usually just the variable/symbol from [symbols.py](symbols.py) and the "lhs" is then the formula
from [checked_functions.py](checked_functions.py).
To display and compute an integral one
uses [``sympy.Integral(..)``](https://docs.sympy.org/latest/modules/integrals/integrals.html) to first display the
desired integral and then ``.doit()`` to actually compute the defined integral.
If you want to compute the integral using an approximation method, there is the option ``"method=quad"`` which needs to
be put into the ``.doit(..)`` call.

## Convert a jupyter notebook to a .py file

Open a console and type ``jupyter nbconvert mynotebook.ipynb --to python``.
