
# Introduction

&nbsp;&nbsp;&nbsp;&nbsp;Thomas Relativistic Electronic Structure Calculation
(TRESC) is used to calculate the electronic structure of non-periodic polyatomic
systems under the Born�COppenheimer approximation, it's designed to do
all-electron single-configuration self-consistent field calculation based on the
static 2-component DKH2 electronic Hamiltonian of a given molecule.<br>
&nbsp;&nbsp;&nbsp;&nbsp;Current program is still under development.

# Algorithms

* Cartesian or spherical-harmonic fragment contracted Gaussian type orbital are
used (input .gbs file).
* Initial guess read from .chk file (by Gaussian HF/KS-SCF) or .ao2mo file
(by TRESC job earlier).
* Symmetric orthogonalisation are used by default, canonical orthogonalisation
will be used to reduce the linear dependence of basis
if the threshold is exceeded.
* Relativistic 1e integrals using RI approach (proposed by Hess *et al.*).
* Non-relativistic 2e integrals using Rys quadrature scheme, screening by
Cauchy-Schwarz scheme.
* Construct Fock matrix via **direct** way, which is time consuming but less
demanding on memory and disk r/w.
* Integrals are consistent with Gaussian program, and relative energies differ
negligibly from Gaussian program.
* If the frontier orbitals are (nearly) degenerate, 2c-SCF results may be far
from initial spin states due to the mixing of alpha
and beta orbitals in degenerate space, try keyword `keepspin` to to avoid this
as possible.
* DIIS(Pulay mixing) can be used to accelerate SCF, dynamic damping can be used
to enhance convergence.
* Basic linear algebra is computed using LAPACK subroutines.
* Both 1e and 2e integrals support OpenMP parallel computation.
* Output $$\left< s^2 \right> \left( L\ddot{o}wdin \right)$$, energy components
and orbital components.
* Support dispersion correction via DFT-D4 program (stand-alone) developed by
Grimme's group.

# Characteristic

## Visualisation of 2-component complex MO

&nbsp;&nbsp;&nbsp;&nbsp;Canonical orbitals will be dumped to a molden format
file containing alpha and beta components with keyword `molden`.
With it, one can generate Gaussian cube format file for any orbital and use
`/scripts/2cvis.py` to visualize it in `mayavi`.<br>
&nbsp;&nbsp;&nbsp;&nbsp;The results are as follows:<br>
<p align="center">
  <img src="docs/figure_2.png" alt="docs/figure_2.png" width="500">
</p>

&nbsp;&nbsp;&nbsp;&nbsp;It's the HOMO of the triplet carbene $$\mathrm{CH}_2$$,
phase deviation from $$\pm {{\mathrm{\pi}}\Bigg/{2}}$$
implies a stronger SOC effect, the plotted results agree with the fact that SOC
intensity is proportional to $${{1}\Bigg/{\mathrm{r}^3}}$$ approximately.

## A special Hamiltonian: SRTP

&nbsp;&nbsp;&nbsp;&nbsp;Second Relativized Thomas Precession (SRTP) is to
conbine the Lorentz vector feature of spin 4-vector
$$\left( 0,\vec{s} \right) $$ and the Lorentz scalar feature of the magnitude of
its spatial components ($$\left| \vec{s} \right|=\hbar /2$$).
'Second Relativized' means the magnitude of spin vector is independent of the
reference frame choice.<br>
&nbsp;&nbsp;&nbsp;&nbsp;To accomplish it, we start with a newly-defined
reference frame transformation rule, which makes the observed $$\vec{s}/s$$
from any frame identical with the observed $$\vec{s}/s$$ from corresponding
frame under the Lorentz transformation rule, but magnitude $$s$$ always
$$\hbar /2$$.<br>
&nbsp;&nbsp;&nbsp;&nbsp;Assuming that frame O' is moving along the x-axis in
frame O, the Lorentz transformation and the newly-defined transformation
lead to different observations.<br>
<p align="center">
  <img src="docs/figure_1.png" alt="docs/figure_1.png" width="500">
</p>

&nbsp;&nbsp;&nbsp;&nbsp;Its mathematical form can be given directly as a
nonlinear equation<br>

$$
\vec{s}\prime=\left( \begin{matrix}
 1&  0&  0&  0\\
 -\gamma \beta _1\zeta&  \left[ 1+\frac{\left( \gamma -1 \right) \beta _{1}^{2}}
 {\beta ^2} \right] \zeta&  \frac{\left( \gamma -1 \right) \beta _1\beta _2}
 {\beta ^2}\zeta&  \frac{\left( \gamma -1 \right) \beta _1\beta _3}{\beta ^2}
 \zeta\\
 -\gamma \beta _2\zeta&  \frac{\left( \gamma -1 \right) \beta _1\beta _2}
 {\beta ^2}\zeta&  \left[ 1+\frac{\left( \gamma -1 \right) \beta _{2}^{2}}
 {\beta ^2} \right] \zeta&  \frac{\left( \gamma -1 \right) \beta _2\beta _3}
 {\beta ^2}\zeta\\
 -\gamma \beta _3\zeta&  \frac{\left( \gamma -1 \right) \beta _1\beta _3}
 {\beta ^2}\zeta&  \frac{\left( \gamma -1 \right) \beta _2\beta _3}{\beta ^2}
 \zeta&  \left[ 1+\frac{\left( \gamma -1 \right) \beta _{3}^{2}}
 {\beta ^2} \right] \zeta\\
\end{matrix} \right) \vec{s}
$$

&nbsp;&nbsp;&nbsp;&nbsp;which $$s_i$$
represent spin components and
$$\beta _i$$
represent velocity components,
$$\gamma$$
represent Lorentz factor
and<br>

$$
\zeta =\sqrt{1+\gamma ^2\beta ^2s^2}
$$

&nbsp;&nbsp;&nbsp;&nbsp;This newly-defined transformation is kinematic, but it
will change the form of Thomas precession dynamically since Thomas precession is
related to the intrinsic property of reference frame transformation.<br>
&nbsp;&nbsp;&nbsp;&nbsp;After some derivation, the contribution of the Thomas
precession to electron energy at low speed can be represented as<br>

$$
H_{\mathrm{SRTP}}=\frac{1}{2}\vec{s}_{\gamma}\cdot
\left( \dot{\vec{\beta}}\times \vec{\beta} \right)
$$

&nbsp;&nbsp;&nbsp;&nbsp;which<br>

$$
s_{\gamma ,i}=\frac{1}{\sqrt{1-\beta _i^{2}}}s_i
$$

&nbsp;&nbsp;&nbsp;&nbsp;Then quantization and use Pauli vector rule to modify
Dirac matrices as<br>

$$
\alpha _i=\left( \frac{\left( 1-\beta _{j}^{2} \right)\left( 1-\beta _{k}^{2}
\right)}{1-\beta _{i}^{2}} \right) ^{\frac{1}{4}}\left( \begin{matrix}
 &  \sigma _i\\
 \sigma _i&  \\
\end{matrix} \right)
$$

&nbsp;&nbsp;&nbsp;&nbsp;This formular leads to the modified electron spinor wave
 function through DKH transformation.<br>
&nbsp;&nbsp;&nbsp;&nbsp;In addition, SRTP effect is of order $$c^{-4}$$, one
have to consider other terms of order $$\geqslant c^{-4}$$ before it, including
radiation effect. Moreover, the lowest order of SRTP still requires the
computation of integrals like $$\langle i|p_{x}^{3}V_{ij}p_y|j\rangle$$, it has
a small effect on results but will significantly increases the one-electron
integral cost.<br>
&nbsp;&nbsp;&nbsp;&nbsp;SRTP is currently has no evidence support, if you're
interested, try keyword `SRTP` when performing DKH2 calculation.<br>

# Upcoming

* support for DFT calculation;
* calculate 2e SOC by SOMF approach;
