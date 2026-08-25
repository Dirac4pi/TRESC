# DKH1-2e Hamiltonian Derivation

For the DKHn Hamiltonian, although the two-electron and one-electron DKH
transformations have the same $c$-order, the former entails much smaller effects
on the relative energy and electronic structure than the latter.
Considering the computational complexity and the fact that higher-order terms
scale as $\mathcal{O}(c^{-4})$, we rigorously truncate the two-electron
Hamiltonian at the $c^{-2}$ order. This means only the free-particle Foldy-Wouthuysen
(fpFW) transformation is applied to the two-electron operators, yielding the
DKH1-2e Coulomb and Breit terms.

$$
\varepsilon_{\mathrm{1,2e},C} = \sum_{i<j}A_{i}A_{j}\left[\frac{1}{r_{ij}} +
R_{i}\frac{1}{r_{ij}}R_{i} + R_{j}\frac{1}{r_{ij}}R_{j}\right]A_{i}A_{j}
$$

$$
\varepsilon_{\mathrm{1,2e},B} = \sum_{i<j}A_{i}A_{j}\left[R_{i}R_{j}\widetilde{B}_{ij}
+ R_{i}\widetilde{B}_{ij}R_{j} + R_{j}\widetilde{B}_{ij}R_{i} +
\widetilde{B}_{ij}R_{i}R_{j}\right]A_{i}A_{j}
$$

where $\widetilde{B}_{ij} =
-\frac{1}{2r_{ij}}\left(\boldsymbol{\sigma}_{i}\cdot\boldsymbol{\sigma}_{j} +
\left(\boldsymbol{\sigma}_{i}\cdot\hat{\boldsymbol{r}}_{ij}\right)\left
(\boldsymbol{\sigma}_{j}\cdot\hat{\boldsymbol{r}}_{ij}\right)\right)$,
$A = \sqrt{\frac{\epsilon+c^{2}}{2\epsilon}}$, and $R =
\frac{c\boldsymbol{\alpha}\cdot\boldsymbol{p}}{\epsilon+c^{2}}$. These equations
are derived from *Transgressing Theory Boundaries: The Generalized Douglas-Kroll
Transformation* by Hess. Dirac's equation indicates that spin-same-orbital
coupling, as the only $c^{-2}$ term in the spin-dependent part of the 2e Coulomb
interaction, occurs only in non-cross terms, thus the first equation lacks the
cross terms between $R_{i}$ and $R_{j}$.

It is noteworthy that for 4-indicator integrals, both equations operate on the
so-called ‘physical representation’ rather than the ‘chemical representation’.
For example:

$$
\left\langle p_{i}V_{ij}p_{j}\right\rangle _{\mu\nu\kappa\lambda} = \left\langle
  \chi_{\mu}(i)\chi_{\nu}(j)\right|\hat{p}_{i}\hat{V}_{ij}\hat{p}_{j}
  \left|\chi_{\kappa}(i)\chi_{\lambda}(j)\right\rangle
$$

and:

$$
\left\langle p_{i}V_{ij}p_{i}\right\rangle _{\mu\nu\kappa\lambda} = \left\langle
  \chi_{\mu}(i)\chi_{\nu}(j)\right|\hat{p}_{i}\hat{V}_{ij}\hat{p}_{i}\left|
    \chi_{\kappa}(i)\chi_{\lambda}(j)\right\rangle
$$

This shows that the operators on both sides operate only on one of the two Gaussian
basis functions on each side, directly yielding 4-indicator integrals that can be
converted into classical and exchange terms of 2-indicator integrals.

Additionally, attention should be paid to how Hess's RI technique can be applied
to the 2-electron DKH integrals. In this case, we need to use the $p^2$ eigenbasis
$\left|\chi^{p2}\right\rangle$ and insert the 2-electron identity:

$$
\begin{aligned}
\left\langle A_{i}A_{j}R_{i}\frac{1}{r_{ij}}R_{i}A_{i}A_{j}\right\rangle
_{\mu\nu\kappa\lambda}^{p2} &= \sum_{mn}\sum_{pq}\left\langle \chi_{\mu}^{p2}(i)
\chi_{\nu}^{p2}(j)\right|A_{i}A_{j}\left|\chi_{m}^{p2}(i)\chi_{n}^{p2}(j)\right\rangle \\
&\quad \times \left\langle \chi_{m}^{p2}(i)\chi_{n}^{p2}(j)\right|R_{i}\frac{1}
{r_{ij}}R_{i}\left|\chi_{p}^{p2}(i)\chi_{q}^{p2}(j)\right\rangle \left\langle
  \chi_{p}^{p2}(i)\chi_{q}^{p2}(j)\right|A_{i}A_{j}\left|\chi_{\kappa}^{p2}(i)
  \chi_{\lambda}^{p2}(j)\right\rangle \\
&= A_{\mu}^{p2}(i)A_{\nu}^{p2}(j)\left\langle \chi_{\mu}^{p2}(i)
\chi_{\nu}^{p2}(j)\right|R_{i}\frac{1}{r_{ij}}R_{i}\left|\chi_{\kappa}^{p2}(i)
\chi_{\lambda}^{p2}(j)\right\rangle A_{\kappa}^{p2}(i)A_{\lambda}^{p2}(j)
\end{aligned}
$$

Since $\left\langle \chi_{\mu}^{p2}(i)\right|A_{i}\left|\chi_{m}^{p2}(i)\right
\rangle = A_{\mu}^{p2}(i)\delta_{\mu m}$ is utilized, the 4-indicator integral
in the $p^2$ eigenbasis in this expression is difficult to handle, so we need to
investigate their contribution to the 2-indicator Coulomb integral:

$$
\begin{aligned}
\left\langle A_{i}A_{j}R_{i}\frac{1}{r_{ij}}R_{i}A_{i}A_{j}\right\rangle
_{\mu\kappa,C}^{p2} &= \sum_{\nu\lambda}\left[\left\langle \chi_{\mu}^{p2}(i)
\chi_{\nu}^{p2}(j)\right|A_{i}A_{j}R_{i}\frac{1}{r_{ij}}R_{i}A_{i}A_{j}\left|
  \chi_{\kappa}^{p2}(i)\chi_{\lambda}^{p2}(j)\right\rangle D_{\nu\lambda}^{p2}\right] \\
&= A_{\mu}^{p2}(i)A_{\kappa}^{p2}(i)\sum_{\nu\lambda}\left[\left\langle
  \chi_{\mu}^{p2}(i)\chi_{\nu}^{p2}(j)\right|R_{i}\frac{1}{r_{ij}}R_{i}\left|
    \chi_{\kappa}^{p2}(i)\chi_{\lambda}^{p2}(j)\right\rangle
    D_{\nu\lambda}^{p2}A_{\nu}^{p2}(j)A_{\lambda}^{p2}(j)\right]
\end{aligned}
$$

Let $D^{p2}$ be the single-electron density matrix in the $p^2$ eigenbasis, and
let $D_{\nu\lambda}^{p2\prime} = D_{\nu\lambda}^{p2}A_{\nu}^{p2}(j)
A_{\lambda}^{p2}(j)$. For the 2-electron integral $X_{ij}$, we have
$\left\langle X_{ij}\right\rangle _{\mu\kappa,C}^{p2} = \left[\Omega^{\dagger}
\left\langle X_{ij}\right\rangle _{C}\Omega\right]_{\mu\kappa}$, so we obtain:

$$
\begin{aligned}
\left\langle A_{i}A_{j}R_{i}\frac{1}{r_{ij}}R_{i}A_{i}A_{j}\right\rangle
_{\mu\kappa,C}^{p2} &= A_{\mu}^{p2}(i)A_{\kappa}^{p2}(i)
\left[\Omega^{\dagger}\left\langle R_{i}\frac{1}{r_{ij}}R_{i}\right\rangle
_{C}^{\prime}\Omega\right]_{\mu\kappa}
\end{aligned}
$$

Similarly, for the 2-indicator exchange integral we get:

$$
\begin{aligned}
\left\langle A_{i}A_{j}R_{i}\frac{1}{r_{ij}}R_{i}A_{i}A_{j}\right\rangle
_{\mu\lambda,X}^{p2} &= A_{\mu}^{p2}(i)A_{\lambda}^{p2}(i)\left[\Omega^{\dagger}
\left\langle R_{i}\frac{1}{r_{ij}}R_{i}\right\rangle _{X}^{\prime}\Omega\right]_{\mu\lambda}
\end{aligned}
$$

To summarize, the complete two-electron Hamiltonian rigorously truncated at the
$c^{-2}$ order (DKH1-2e) is:

$$
H_{\mathrm{DKH1,2e}} = \varepsilon_{\mathrm{1,2e},C} + \varepsilon_{\mathrm{1,2e},B}
$$

# Contribution of 4-Index Integrals to the Dirac-Coulomb Matrix

Kinematic operators:
$$A_p = \sqrt{\frac{E_p + m c^2}{2 E_p}}, \quad R_p = \frac{c \, \mathbf{\alpha}
\cdot \mathbf{p}}{E_p + m c^2} = \mathbf{\alpha} \cdot \mathbf{P}_p = \mathcal{R}_p
\, \mathbf{\alpha} \cdot \mathbf{p}$$

Two-electron Coulomb Hamiltonian ($H_1^{2e}(V_C)$), subscripts distinguish between
electron 1 and electron 2:
$$H_1^{2e}(V_C) = \sum_{i < j} A_i A_j \left\{ \frac{e^2}{r_{ij}} + R_i
\frac{e^2}{r_{ij}} R_i + R_j \frac{e^2}{r_{ij}} R_j + R_i R_j \frac{e^2}{r_{ij}}
R_i R_j \right\} A_i A_j$$
$$= (A_p)_1 (A_p)_2 \left[ \frac{1}{r_{12}} + (\mathcal{R}_p)_1 (\mathbf{\alpha}_1
\cdot \mathbf{p}_1) \frac{1}{r_{12}} (\mathbf{\alpha}_1 \cdot \mathbf{p}_1)
(\mathcal{R}_p)_1 + (\mathcal{R}_p)_2 (\mathbf{\alpha}_2 \cdot \mathbf{p}_2)
\frac{1}{r_{12}} (\mathbf{\alpha}_2 \cdot \mathbf{p}_2) (\mathcal{R}_p)_2 \right]
(A_p)_1 (A_p)_2$$
$$= (A_p)_1 (A_p)_2 \frac{1}{r_{12}} (A_p)_1 (A_p)_2 + (A_p \mathcal{R}_p)_1
(A_p)_2 (\mathbf{\alpha}_1 \cdot \mathbf{p}_1) \frac{1}{r_{12}} (\mathbf{\alpha}_1
\cdot \mathbf{p}_1) (A_p \mathcal{R}_p)_1 (A_p)_2 + (A_p)_1 (A_p \mathcal{R}_p)_2
(\mathbf{\alpha}_2 \cdot \mathbf{p}_2) \frac{1}{r_{12}} (\mathbf{\alpha}_2 \cdot
\mathbf{p}_2) (A_p)_1 (A_p \mathcal{R}_p)_2$$

Term-by-term expansion for integral $\langle i k | j l \rangle$, superscripts
distinguish electron 1 and electron 2:
* **Term 1**: $(A_p)_i^1 (A_p)_k^2 \frac{1}{r_{12}} (A_p)_j^1 (A_p)_l^2$
* **Term 2**: $(A_p \mathcal{R}_p)_i^1 (A_p)_k^2 (\mathbf{\sigma}_1 \cdot
\mathbf{p}_1) \frac{1}{r_{12}} (\mathbf{\sigma}_1 \cdot \mathbf{p}_1) (A_p
\mathcal{R}_p)_j^1 (A_p)_l^2$
* **Term 3**: $(A_p)_i^1 (A_p \mathcal{R}_p)_k^2 (\mathbf{\sigma}_2 \cdot
\mathbf{p}_2) \frac{1}{r_{12}} (\mathbf{\sigma}_2 \cdot \mathbf{p}_2) (A_p)_j^1
(A_p \mathcal{R}_p)_l^2$

Composite Density Matrices, superscripts distinguish electron 1 and electron 2:
$$D_1^1 = (A_p)_k D (A_p)_l \text{for term 1 and term 2}$$
$$D_2^1 = (A_p R_p)_k^2 D (A_p R_p)_l^2 \text{for term 3}$$

Expanded Pauli Operator Product (unsimplified form):
$$(\sigma_x P_x + \sigma_y P_y + \sigma_z P_z) V (\sigma_x P_x + \sigma_y P_y +
\sigma_z P_z)$$
$$= P_x V P_x \sigma_x^2 + P_x V P_y \sigma_x \sigma_y + P_x V P_z \sigma_x
\sigma_z + P_y V P_x \sigma_y \sigma_x + P_y V P_y \sigma_y^2 + P_y V P_z
\sigma_y \sigma_z + P_z V P_x \sigma_z \sigma_x + P_z V P_y \sigma_z \sigma_y +
P_z V P_z \sigma_z^2$$


## Four-Index Integral Contributions, Spin Routing, and Symmetry Deduplication

Using the chemists' notation, we define the scalar and general two-component
momentum integrals as:
$$
(ij|kl) \equiv \langle \chi_i(\mathbf{r}_1)\chi_j(\mathbf{r}_2) |r_{12}^{-1}|
\chi_k(\mathbf{r}_1)\chi_l(\mathbf{r}_2)\rangle
$$
$$
X_{ijkl}^{ab} \equiv (P_a i,j|P_b k,l), \qquad a,b\in\{x,y,z\}
$$
Owing to the permutation symmetry of the two-electron interaction, a single
computed four-index integral contributes to multiple Coulomb ($J$) and Exchange
($K$) matrix elements.

---

### 1. Coulomb Contraction & Spin Routing

The general Coulomb contraction traces over the spectator electron (electron 2):
$$
J_{ij}^{\sigma_1\sigma_3} = \sum_{kl} \sum_{\sigma_2\sigma_4}
X_{ijkl}^{\sigma_1\sigma_2;\sigma_3\sigma_4} D_{lk}^{\sigma_4\sigma_2}
$$

Depending on the operator location in the $c^{-2}$ DKH1 terms
($V_C$, $R_1V_CR_1$, $R_2V_CR_2$), the spin routing behaves differently:

*   **Scalar ($V_C$) and Operator on Electron 1 ($R_1V_CR_1$):**
    The spectator electron contracts with the total scalar density
    $D^T_{lk} = D_{lk}^{\alpha\alpha} + D_{lk}^{\beta\beta}$. For $R_1V_CR_1$,
    the Pauli matrices dictate the spin destination. The physical spin-dependent
    (SD) antisymmetric combination is:
    $$
    i\left[ (P_x i,j|P_y k,l) - (P_y i,j|P_x k,l) \right]
    $$
*   **Operator on Electron 2 ($R_2V_CR_2$):**
    The Pauli matrices act on the traced electron, absorbing the operator into a
    spin-dependent density $\Gamma$:
    $$
    \Gamma_{lk}^{ab} = \sum_{\rho\tau} (\sigma_a\sigma_b)_{\rho\tau} D_{lk}^{\tau\rho}
    $$
    The resulting Coulomb matrix acts as a pure scalar to electron 1:
    $$
    J_{ij}^{\sigma\sigma'} = \delta_{\sigma\sigma'} \sum_{kl,ab} (P_a i,j|P_b k,l)
    \Gamma_{lk}^{ab}
    $$

---

### 2. Exchange Contraction

The exchange contraction connects both electron lines, resulting in a distinct
spin routing:
$$
K_{ik}^{\sigma\tau} = \sum_{jl} \sum_{\rho\lambda}
X_{ijkl}^{\sigma\rho;\lambda\tau} D_{lj}^{\lambda\rho}
$$
Unlike the Coulomb matrix, the exchange matrix does not require $\sigma=\tau$.
Off-diagonal density blocks ($D^{\alpha\beta}$, $D^{\beta\alpha}$) actively
contribute depending on the Pauli operators bridging the exchanged indices.

---

### 3. Permutation Symmetry and Deduplication

A single scalar integral has four canonical Fock-matrix destinations:
$$
(ij|kl) \Rightarrow \{J_{ij}, J_{kl}, K_{ik}, K_{jl}\}
$$

For spin-dependent integrals $X_{ijkl}^{ab} = (P_a i,j|P_b k,l)$, the pair
interchange symmetry $(P_a i,j|P_b k,l) \leftrightarrow (P_b k,l|P_a i,j)$ also
transforms the Cartesian labels:
*   **$a \neq b$ (Cross Components):** $(P_a i,j|P_b k,l)$ and $(P_b k,l|P_a i,j)$
are distinct tensor components. Retain the ordered Cartesian pair.
*   **$a = b$ (Co-directional):** The integral retains standard pair-interchange
symmetry.

**Core Deduplication Rule:**
A factor of $1/2$ is required *only* when two nominally distinct Fock-matrix
destinations generated from the same integral collapse onto the exact same matrix
element. **Deduplicate equivalent Fock destinations, not simply equal AO indices.**

$$
J_{ij} = J_{kl} \iff (i,j) = (k,l)
$$
$$
K_{ik} = K_{jl} \iff (i,k) = (j,l)
$$

---

### 4. Summary of the Complete Routing Logic

The routing rules for the $c^{-2}$ DKH1 two-electron Coulomb terms are summarized
as follows:

**Coulomb Destinations:**
$$
\begin{array}{c|c|c}
\text{Operator Location} & \text{Contracted Density} & \text{Spin Destination} \\
\hline
V_C & D^T & J \propto I \\
R_1V_CR_1 & D^T & J \propto \sigma_a\sigma_b \\
R_2V_CR_2 & \Gamma^{ab} & J \propto I
\end{array}
$$

**Exchange Destinations:**
$$
\begin{aligned}
V_C &: \quad D^{\tau\sigma} \rightarrow K^{\sigma\tau} \\
R_1V_CR_1 &: \quad D^T \rightarrow (\sigma_a\sigma_b)\,K \\
R_2V_CR_2 &: \quad D^{\sigma\rho} \text{ contracts with } (\sigma_a\sigma_b)_{\rho\tau}
\end{aligned}
$$

**Algorithmic Workflow:**
$$
\text{Integral Symmetry} \longrightarrow \text{Fock Destinations} \longrightarrow
\text{Spin Routing} \longrightarrow \text{Deduplication}
$$