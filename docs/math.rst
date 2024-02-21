Definitions 
===========
Let $x_{j} = (x_{j,1}, x_{j,2}) \in \mathbb{R}^{2}$, $j=1,2,\ldots N$, denote a collection
of source locations and let $x = (x_{1}, x_{2}) \in \mathbb{R}^{2}$ denote a
target location. Unless stated otherwise, the gradients are with
respect to the target variable $x$.

Laplace FMM
*************
The Laplace FMM comes in three varieties

* rfmm2d: charges, dipole strengths, potentials, their
  gradients, and Hessians are all double precision 
* lfmm2d: all of the above quantities are double complex
* cfmm2d: all quantities are double complex, except there are no dipole orientation vectors. Instead, the dipole kernel is assumed to be of the form $1/(z - \xi_{j})$ where $z, \xi_{j}$ are the target and source locations, viewed as points in the complex
  plane given by $z = x_{1} + i\cdot x_{2}$, and $\xi_{j} = x_{j,1} + i \cdot x_{j,2}$. 
  In this case, instead of gradients and Hessians, the first and second derivatives of 
  the potential with respect to $z$ are returned and, with a slight abuse of
  notation, we set $\frac{\mathrm{d}}{\mathrm{d} z} \log(\|z\|) = 1/z$.
 

Laplace FMM (rfmm2d)
---------------------
Let $c_{j} \in \mathbb{R}$, 
$j=1,2,\ldots N$, 
denote a collection of charge strengths, let $v_{j} \in \mathbb{R}$,
$j=1,2,\ldots N$, 
denote a collection of dipole strengths, and let $d_{j} \in \mathbb{R}^{2}$,
$j=1,2,\ldots N$, denote the corresponding dipole orientation vectors.

The Laplace FMM (rfmm2d) computes 
the potential $u(x) \in \mathbb{R}$
given by

.. math::
   :label: rlap_nbody

    u(x) = \sum_{j=1}^{N} c_{j} \log{(\|x-x_{j}\|)} - v_{j} d_{j}\cdot \nabla \log{(\|x-x_{j}\|)}  \, , 

and its gradient $\nabla u(x) \in \mathbb{R}^{2}$
at the source and target locations. When $x=x_{j}$, the term
corresponding to $x_{j}$ (the self-interaction) is omitted from the sum.


Laplace FMM (lfmm2d)
---------------------
Let $c_{j} \in \mathbb{C}$, 
$j=1,2,\ldots N$, 
denote a collection of charge strengths, let $v_{j} \in \mathbb{C}$,
$j=1,2,\ldots N$, 
denote a collection of dipole strengths, and let $d_{j} \in \mathbb{R}^{2}$,
$j=1,2,\ldots N$, denote the corresponding dipole orientation vectors.

The Laplace FMM (rfmm2d) computes 
the potential $u(x) \in \mathbb{C}$ 
given by

.. math::
   :label: llap_nbody

    u(x) = \sum_{j=1}^{N} c_{j} \log{(\|x-x_{j}\|)} - v_{j} d_{j}\cdot \nabla \log{(\|x-x_{j}\|)}  \, , 

and its gradient $\nabla u(x) \in \mathbb{C}^{2}$
at the source and target locations. (The outputs are complex-valued in this case because
the $\{ c_j \}$ and $\{ v_j \}$ are complex.)
When $x=x_{j}$, the term
corresponding to $x_{j}$ is omitted from the sum.


Laplace FMM (cfmm2d)
---------------------
Let $c_{j} \in \mathbb{C}$, 
$j=1,2,\ldots N$, 
denote a collection of charge strengths, and let $v_{j} \in \mathbb{C}$,
$j=1,2,\ldots N$, 
denote a collection of dipole strengths.

For $z, \xi_{j} \in \mathbb{C}$, 
the Laplace FMM (cfmm2d) computes 
the potential $u(z) \in \mathbb{C}$ 
given by

.. math::
   :label: clap_nbody

    u(z) = \sum_{j=1}^{N} c_{j} \log{(\|z-\xi_{j}\|)} - \frac{v_{j}}{z-\xi_{j}}  \, , 

and its 
derivatives $u'(z), u''(z) \in
\mathbb{C}$
at the source and target locations. When $z=\xi_{j}$, the term
corresponding to $\xi_{j}$ is omitted from the sum. 
As noted above, we define $\frac{\mathrm{d}}{\mathrm{d} z} \log(\|z\|) = 1/z$.


Helmholtz FMM
*************
Let $c_{j} \in \mathbb{C}$, 
$j=1,2,\ldots N$, 
denote a collection of charge strengths, let $v_{j} \in \mathbb{C}$,
$j=1,2,\ldots N$, 
denote a collection of dipole strengths, and let $d_{j} \in \mathbb{R}^{2}$,
$j=1,2,\ldots N$, denote the corresponding dipole orientation vectors.
Let $k\in\mathbb{C}$ denote the wave number (the Helmholtz 
parameter). 

The Helmholtz FMM computes 
the potential $u(x) \in \mathbb{C}$ 
given by

.. math::
   :label: helm_nbody

    u(x) = \sum_{j=1}^{N} c_{j} H_{0}^{(1)}(k\|x-x_{j}\|) - v_{j} d_{j}\cdot \nabla H_{0}^{(1)}(k\|x-x_{j}\|)  \, , 

and its gradient $\nabla u(x) \in \mathbb{C}^{2}$
at the source and target locations, where $H_{0}^{(1)}$ is the Hankel function
of the first kind of order $0$. When $x=x_{j}$, the term
corresponding to $x_{j}$ is omitted from the sum.


Biharmonic FMM
***************
Let $c_{j} = (c_{j,1}, c_{j,2})\in \mathbb{C}^2$, 
$j=1,2,\ldots N$, 
denote a collection of charge strengths, and let
$v_{j} = (v_{j,1}, v_{j,2}, v_{j,3}) \in \mathbb{C}^{3}$,
$j=1,2,\ldots N$, 
denote a collection of dipole strengths.

For $z, \xi_j \in \mathbb{C}$, the biharmonic FMM computes 
the potential $u(z)$ and its `gradient` = 
$(P_{z} \frac{\mathrm{d}}{\mathrm{d}z}, P_{\overline{z}} \frac{\mathrm{d}}{\mathrm{d}z}, \frac{\mathrm{d}}{\mathrm{d}\overline{z}})$
given by

.. math::
   :label: biharm_nbody

    u(z) &= \sum_{j=1}^{N} 2 \, c_{j,1} \log{\|z - \xi_{j}\|} + 
    c_{j,2} \frac{z-\xi_{j}}{\overline{z-\xi_{j}}} + \frac{v_{j,1}}{z - \xi_{j}} + 
    \frac{v_{j,3}}{\overline{z-\xi_{j}}} + 
    v_{j,2} \frac{z - \xi_{j}}{(\overline{z-\xi_{j}})^2} \, , \\
    P_{z} \frac{\mathrm{d}}{\mathrm{d} z}u(z) &= \sum_{j=1}^{N} \frac{c_{j,1}}{z-\xi_{j}} -
    \frac{v_{j,1}}{(z-\xi_{j})^2} \, \\
    P_{\overline{z}} \frac{\mathrm{d}}{\mathrm{d} z} u(z) &= 
    \sum_{j=1}^{N} \frac{c_{j,2}}{\overline{z-\xi_{j}}} + 
    \frac{v_{j,2}}{(\overline{z-\xi_{j}})^2} \,
    ,\\
    \frac{\mathrm{d}}{\mathrm{d}\overline{z}} u(z) &= 
    \sum_{j=1}^{N} \frac{c_{j,1}}{\overline{z-\xi_{j}}} - 
    c_{j,2} \frac{z-\xi_{j}}{(\overline{z-\xi_{j}})^2} - 
    \frac{v_{j,3}}{(\overline{z-\xi_{j}})^2} - 
    2v_{j,2} \frac{z - \xi_{j}}{(\overline{z-\xi_{j}})^3} \, , \\

at the source and target locations. 
When $z=\xi_{j}$, the term
corresponding to $\xi_{j}$ is omitted from the sum. 
The expression $P_{z} \frac{\mathrm{d}}{\mathrm{d}z}$
denotes the component of the derivative $\frac{\mathrm{d}}{\mathrm{d} z}$
which is purely a function of $z$, and the expression
$P_{\overline{z}} \frac{\mathrm{d}}{\mathrm{d}z}$
denotes the component of the derivative $\frac{\mathrm{d}}{\mathrm{d} z}$
which is purely a function of $\overline{z}$.


Modified Biharmonic FMM
************************

Let $G^{\textrm{mbh}}(x,y)$ denote the modified biharmonic
Green's function given by

.. math::
    G^{\textrm{mbh}}(x,y) = \frac{1}{2\pi \beta^2}\left(K_{0}(\beta \|x-y \|) - \log{\|x-y\|}\right)

where $K_{0}$ is the modified Bessel function of order $0$, and $\beta$ is the
modified biharmonic wavenumber.


Let $c_{j} \in \mathbb{R}$, 
denote a collection of charge strengths, let
$v_{j} \in \mathbb{R}$, 
denote a collection of dipole strengths, let
$d_{j} = (d_{j,1}, d_{j,2})$ denote a collection
of dipole vectors, let
$q_{j} \in \mathbb{R}$ denote a collection of 
quadrupole strengths, let
$w_{j} = (w_{j,1}, w_{j,2}, w_{j,3}) \in \mathbb{R}^{3}$, 
denote a collection of quadrupole three-vectors, let
$o_{j}$ denote a collection of octopole strengths, and let
$p_{j} = (p_{j,1}, p_{j,2}, p_{j,3}, p_{j,4}) \in \mathbb{R}^{4}$, 
denote a collection of octopole four-vectors.

The modified biharmonic FMM computes the potential $u(x)\in \mathbb{R}$ 
given by

.. math::
   :label: modbiharm_nbody

    u(x) = \sum_{j=1}^{N} &c_{j} G^{\textrm{mbh}}(x,x_{j}) + 
    v_{j} d_{j} \cdot \nabla_{y} G^{\textrm{mbh}}(x,x_{j}) + \\
    &q_{j} \left(w_{j,1} \partial_{y_{1},y_{1}} G^{\textrm{mbh}}(x,x_{j}) + w_{j,2}
    \partial_{y_{1},y_{2}} G^{\textrm{mbh}}(x,x_{j}) + w_{j,3}
    \partial_{y_{2},y_{2}} G^{\textrm{mbh}}(x,x_{j}) \right) + \\
    &o_{j} \big( p_{j,1} \partial_{y_{1},y_{1},y_{1}} G^{\textrm{mbh}}(x,x_{j}) +
    p_{j,2} \partial_{y_{1},y_{1},y_{2}} G^{\textrm{mbh}}(x,x_{j}) + \\
    &p_{j,3} \partial_{y_{1},y_{2},y_{2}} G^{\textrm{mbh}}(x,x_{j}) +
    p_{j,4} \partial_{y_{2},y_{2},y_{2}} G^{\textrm{mbh}}(x,x_{j}) \big) \, ,

and its gradients $\nabla u(x) \in \mathbb{R}^{2}$
at the source and target locations. When $x=x_{j}$, the term
corresponding to $x_{j}$ is omitted from the sum.

Stokes FMM
************

Let $G^{\textrm{stok}}(x,y)$ denote the Stokeslet given by

.. math::
    G^{\textrm{stok}}(x,y) = \frac{1}{2}
    \begin{bmatrix}
    -\log{\|x-y \|} +  \frac{(x_{1}-y_{1})^2}{\|x-y\|^2} & \frac{(x_{1}-y_{1})
    (x_{2}-y_{2})}{\|x-y \|^2} \\
    \frac{(x_{1}-y_{1})(x_{2}-y_{2})}{\|x-y \|^2} &
    -\log{\|x-y \|} +  \frac{(x_{2}-y_{2})^2}{\|x-y \|^2} 
    \end{bmatrix} \, ,

and let $P^{\textrm{stok}}(x,y)$ denote the associated pressure tensor

.. math::
    P^{\textrm{stok}}(x,y) = \frac{1}{\|x-y \|^2}\begin{bmatrix}
    (x_{1}-y_{1}) &
    (x_{2}-y_{2})
    \end{bmatrix} \, .

Let $T^{\textrm{stok}}(x,y)$ denote the Stresslet whose action on a vector $v$
is given by

.. math::
    v\cdot T^{\textrm{stok}}(x,y) = -\frac{2 v \cdot (x-y)}{\|x-y \|^4}
    \begin{bmatrix}
    (x_{1} - y_{1})^2 & (x_{1}-y_{1})(x_{2}-y_{2}) \\
    (x_{1}-y_{1})(x_{2}-y_{2}) & (x_{2} - y_{2})^2
    \end{bmatrix} \, ,

and let $\Pi^{\textrm{stok}} (x,y)$ denote its associated pressure tensor given by

.. math::
    v\cdot \Pi(x,y)^{\textrm{stok}} = -\frac{v}{\|x-y \|^2} + \frac{2 v \cdot(x-y)}{\|x-y \|^4}
    \begin{bmatrix}
    (x_{1}-y_{1}) \\
    (x_{2}-y_{2})
    \end{bmatrix} \, .

Let $c_{j} \in \mathbb{R}^2$, 
$j=1,2,\ldots N$, 
denote a collection of Stokeslet strengths, let $v_{j} \in \mathbb{R}^2$,
$j=1,2,\ldots N$, 
denote a collection of Stresslet strengths, and let $d_{j} \in \mathbb{R}^{2}$,
$j=1,2,\ldots N$, denote the corresponding Stresslet orientation vectors.

The Stokes FMM computes 
the potential $u(x) \in \mathbb{R}^2$, its gradient $\nabla u(x) \in
\mathbb{R}^{2\times 2}$, and the pressure $p$ given by

.. math::
   :label: stok_nbody

    u(x) &= \sum_{j=1}^{N} G^{\textrm{stok}}(x,x_{j}) c_{j} + d_{j} \cdot
    T^{\textrm{stok}}(x,x_{j}) \cdot v_{j} \, , \\
    p(x) &= \sum_{j=1}^{N} c_{j} P^{\textrm{stok}}(x,x_{j}) + d_{j} \cdot
    \Pi^{\textrm{stok}}(x,x_{j}) \cdot v_{j}^{T}

at the source and target locations. When $x=x_{j}$, the term
corresponding to $x_{j}$ is omitted from the sum.

Vectorized versions   
*******************
The vectorized versions of the FMMs compute the same sums
as above for a set of problems in which the source and target
locations are constant but multiple values of the charges, dipoles, etc 
are specified at each source. Given a set of problems with this structure,
the vectorized versions are faster than calling the standard
FMM multiple times in sequence. 

For example, let $c_{\ell,j}\in\mathbb{C}$, 
$j=1,2,\ldots N$, $\ell=1,2,\ldots n_{d}$
denote a collection of $n_{d}$ charge strengths, and
let $v_{\ell,j} \in \mathbb{C}$, $d_{\ell,j} \in \mathbb{R}^2$ 
denote a collection of $n_{d}$ dipole strengths and orientation vectors. 
Then the vectorized Helmholtz FMM computes the potentials $u_{\ell}(x) \in
\mathbb{C}$
given by

.. math::
    :label: helm_nbody_vec

    u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} H_{0}^{(1)}(k\|x-x_{j}\|) - v_{\ell,j} d_{\ell,j}\cdot \nabla H_{0}^{(1)}(k\|x-x_{j}\|)  \, , 

and its gradients $\nabla u_{\ell}(x) \in \mathbb{C}^{2}$
at the source and target locations. 

.. note::

   In double precision arithmetic, two numbers which are
   within machine precision of each other cannot be
   distinguished. In order to account for this, suppose that the sources
   and targets are contained in a square with side length $L$, then
   for all $x$ such that $\| x-x_{j} \| \leq L \varepsilon_{\textrm{mach}}$,
   the term corresponding to $x_{j}$ is omitted from the sum.
   Here $\varepsilon_{\textrm{mach}} = 2^{-52}$ is machine precision.

