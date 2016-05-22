NormalForm
==========

| Explicitly computes the smooth transformation that maps a complicated dynamical system (near a local bifurcation) to a simpler low-dimensional normal form system that preserves the same local behavior.

Overview
--------

This is a package for Mathematica. 

To install, unpack into the location given by ``$UserBaseDirectory``.

| For more information on these algorithms for deterministic systems:
| J. Murdock (2003) "Normal Forms and Unfoldings for Local Dynamical Systems"
|
| For application to systems of stochastic differential equations:
| Aburn, Holmes, Daffertshofer and Breakspear "Normal form transformations explain effect of noise near Hopf bifurcations" (in prep.)

functions
---------

| ``NormalFormTransformation[rhs, {x1,...,xn}, {u1,...,un}, m]`` transforms the dynamical system with right hand side ``rhs`` (expressed in original variables ``{xi}``) to a simpler system (normal form to order ``m``) in the new variables ``{ui}``. Returns a pair ``{newrhs, trans}`` where ``newrhs`` is the transformed system and ``trans`` is a smooth invertible coordinate transformation that maps ``rhs`` to ``newrhs``. N.B. it is assumed that the linear part of the system has already been transformed to Jordan real form.
| Options:
| ``Verbose->True`` will cause it to print out working at each step.
| ``BifurcationParameters->{eps1, eps2, ...}`` set which symbols in ``rhs`` should be interpreted as bifurcation parameters. Default is ``Global`\[Epsilon]``.
| ``AsymptoticScaling->{symbol1^exponent1,...}`` advise what asymptotic scaling to assume when truncating the resulting power series. The default is ``{x1,...,xn,Sqrt[\[Epsilon]]}`` (which means that \[Epsilon] is taken to be the same order as the x_i squared).
| ``Augmented``: whether to compute the normal form of the augmented system, that is with phase space extended with dimensions for the (rescaled) bifurcation parameters and their equations \dot{\alpha}=0, deriving a transformation dependent on the bifurcation parameters. If False (the default) then the normal form will be found with respect to the dynamical variables only.
|
|
| ``TransformNoisyHopf[rhs, {x1,...,xn}, {\[Sigma]1,...,\[Sigma]n}, {\[Xi]1,...\[Xi]n}, r, {new\[Xi]1, new\[Xi]2}]`` takes the stochastic dynamical system with right hand side ``rhs`` (expressed in variables ``{xi}``, small noise parameters ``{\[Sigma]i}`` and Langevin noise symbols ``{\[Xi]i}`` with Stratonovich interpretation of any multiplicative noise) and transforms it to a simple circular 2 dimensional Hopf normal form system (expressed in new polar variables ``{r, \[Theta]}`` and new Langevin noise symbols ``{new\[Xi]1, new\[Xi]2}``). N.B. It is assumed that the linear part of the system has already been transformed to Jordan real form, with Hopf bifurcation in first two variables at the origin.
| Options: 
| ``Verbose->True`` (as above)
| ``BifurcationParameters->{eps}`` (as above)
| ``AsymptoticScaling->{symbol1^exponent1,...}`` with default value ``{x1,...,xn,Sqrt[\[Epsilon],\[Sigma]1,...\[Sigma]n}`` (i.e. by default the noise strengths \[Sigma] are taken to be of the same order as the x_i when truncating the resulting power series)
| ``MaxOrder->n`` The default value 3 (compute all terms, including noise effects, up to third order) is sensible, unless the Hopf is degenerate.
| ``Augmented``: whether to compute the normal form of the augmented system, that is with phase space extended with dimensions for the (rescaled) bifurcation parameters and their equations \dot{\alpha}=0, deriving a transformation dependent on the bifurcation parameters. If False (the default) then the normal form will be found with respect to the dynamical variables only.
| ``Average``: whether to average around the cycle. (default True)
| ``Rescale``: whether to linearly rescale the radial variable to make the coefficient of the R^3 term -1. (default True)
|
|
| ``MultiSeries[vectorField, {x1^exp1, ...}, maxOrder]`` Multivariate power series, supporting different asymptotic scaling of variables.
|
|
| ``TransformContravariant[U, R]`` applies the near-identity coordinate transformation ``U`` to transform the contravariant vector field ``R(u)``. Both ``U`` and ``R`` should be given in the form of ``MultiSeries``.
|
|
| ``BalanceMatrix[A]`` returns the pair ``{T, B}`` where ``T`` is a similarity transformation and ``B`` is the transformed matrix, ``B = T^-1.A.T``, such that ``B`` is as close to symmetric as possible. This is used to improve an ill-conditioned matrix ``A``, allowing eigenvalues and eigenvectors to be computed more precisely from matrix ``B``. Ref: Parlett and Reinsch (1969)

TODO
----
- Should also automate the preparation step (translation and linear transformation to put the linear part of system in Jordan real form with Hopf at the origin in the first two variables). Currently you need to do this preparation separately before using ``NormalFormTranformation[]`` or ``TransformNoisyHopf[]``.

- Implement the more efficient Lie algebra based normal form algorithms given in Murdock (2003).
