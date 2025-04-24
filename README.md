# Urfit

Simple stats package using some things from GSL (but not the optimisation, that is home-spun). The code will do either bootstrap or (bias-correcting double) jackknife on various files my other codes output. It then will perform (fully-Correlated or Uncorrelated) simultaneous fits if the user desires. It computes effective masses of correlators and all sorts of fancy time-series things depending on the analysis required. Also supports the addition of Bayesian priors for fit parameters.

## Minimisers

The code has various minimisers built-in for the \chi^2 functional. The best and fastest for most jobs will be the Marquardt-Levenburg

    LM

And there is also a Conjugate Gradient

    CG

method for those who want to try something worse. Never use the minimisation SD, although it exists.

If the function we have is too complicated for the specification of derivatives I have the following methods

    POWELL, SIMPLEX, GA

Implementing Powell's method from Numerical Recipes, the downhill simplex (most closely resembles SciPy's implementation), or a genetic algorithm with tournament selection and elitism which is just fun to play around with. These only require function evaluations and no derivatives.

If we are lucky enough to have a problem that is linear in the parameters, i.e. can be expressed as

   y = f[0] + f[1]*x^n + f[2]*x^m

With n and m being real numbers we can skip minimisation altogether and directly solve this Generalised Least Squares problem using the routine

     GLS

Provided the user has specified the matrix form needed to invert. Take a look at the polynomial fit routine for an indication of what I mean.

All of these minimisers will perform fits in parallel over bootstraps and jackknifes (except Powell's because the NR interface I followed is miserable).

## Simultaneous by design

The idea of this code was to provide natural sharing of simultaneous fit parameters. Each fit logically orders their parameters, e.g. a single exponential fit has two

    y = p[0]*exp(-p[1]*x)

If we set p[1], the exponent, to be simultaneous over two data sets then the fit will have 3 fit parameters and be ordered

   y_1 = p[0]*exp(-p[1]*x)
   y_2 = p[2]*exp(-p[1]*x)

This form is true for most of the basic fit functions specified, some of the more nuanced and complicated versions sometimes change this and for that you are on your own.