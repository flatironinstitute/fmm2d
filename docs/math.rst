Definitions 
===========
Let $x_{j} = (x_{j,1}, x_{j,2}) \in \mathbb{R}^{2}$, $i=1,2,\ldots N$, denote a collection
of source locations and let $x = (x_{1}, x_{2}) \in \mathbb{R}^{2}$ denote a
target location. Unless stated otherwise, the gradients are with
respect to the target variable $x$.

Laplace FMM
*************
The Laplace FMM comes in three varieties

* rfmm2d: Charges, dipole strengths, potentials, their
  gradients, and Hessians are all double precision 
* lfmm2d: All of the above quantities are double complex
* cfmm2d: Same as lfmm2d except that there are no dipole orientation vectors,
  and the dipole kernel is given by $1/(z - \xi_{j})$ where $z, \xi_{j}$
  are the target, and source locations viewed as points in the complex
  plane given by
  $z = x_{1} + i\cdot x_{2}$, and $\xi_{j} = x_{j,1} + i \cdot x_{j,2}$. 
  Moreover, the gradients, and Hessians
  in this case are derivatives with respect to $z$, and in a slight abuse of
  notation, we set $\frac{\mathrm{d}}{\mathrm{d} z} \log(\|z\|) = 1/z$
 

Laplace FMM (rfmm2d)
---------------------
Let $c_{j} \in \mathbb{R}$, 
$j=1,2,\ldots N$, 
denote a collection of charge strengths, $v_{j} \in \mathbb{R}$,
$j=1,2,\ldots N$, 
denote a collection of dipole strengths, and $d_{j} \in \mathbb{R}^{2}$,
$j=1,2,\ldots N$, denote the corresponding dipole orientation vectors.

The Laplace FMM (rfmm2d) computes 
the potential $u(x) \in \mathbb{R}$
given by

.. math::
   :label: rlap_nbody

    u(x) = \sum_{j=1}^{N} c_{j} \log{(\|x-x_{j}\|)} - v_{j} d_{j}\cdot \nabla \log{(\|x-x_{j}\|)}  \, , 

and its gradient $\nabla u(x) \in \mathbb{R}^{2}$
at the source and target locations. When $x=x_{j}$, the term
corresponding to $x_{j}$ is dropped from the sum.


Laplace FMM (lfmm2d)
---------------------
Let $c_{j} \in \mathbb{C}$, 
$j=1,2,\ldots N$, 
denote a collection of charge strengths, $v_{j} \in \mathbb{C}$,
$j=1,2,\ldots N$, 
denote a collection of dipole strengths, and $d_{j} \in \mathbb{R}^{2}$,
$j=1,2,\ldots N$, denote the corresponding dipole orientation vectors.

The Laplace FMM (rfmm2d) computes 
the potential $u(x) \in \mathbb{C}$ 
given by

.. math::
   :label: llap_nbody

    u(x) = \sum_{j=1}^{N} c_{j} \log{(\|x-x_{j}\|)} - v_{j} d_{j}\cdot \nabla \log{(\|x-x_{j}\|)}  \, , 

and its gradient $\nabla u(x) \in \mathbb{C}^{2}$
at the source and target locations. When $x=x_{j}$, the term
corresponding to $x_{j}$ is dropped from the sum.


Laplace FMM (cfmm2d)
---------------------
Let $c_{j} \in \mathbb{C}$, 
$j=1,2,\ldots N$, 
denote a collection of charge strengths, and $v_{j} \in \mathbb{C}$,
$j=1,2,\ldots N$, 
denote a collection of dipole strengths.

The Laplace FMM (cfmm2d) computes 
the potential $u(z) \in \mathbb{C}$ 
given by

.. math::
   :label: clap_nbody

    u(z) = \sum_{j=1}^{N} c_{j} \log{(\|z-\xi_{j}\|)} - \frac{v_{j}}{z-\xi_{j}}  \, , 

and its gradient $\frac{\mathrm{d}}{\mathrm{d}z} u(z) \in
\mathbb{C}$
at the source and target locations. When $z=\xi_{j}$, the term
corresponding to $\xi_{j}$ is dropped from the sum. Here
$z = x(1) + i x(2)$, and $\xi_{j} = x_{j,1} + ix_{j,2}$, are the complex
numbers corresponding to $x$ and $x_{j}$ respectively.


Helmholtz FMM
*************
Let $c_{j} \in \mathbb{C}$, 
$j=1,2,\ldots N$, 
denote a collection of charge strengths, $v_{j} \in \mathbb{C}$,
$j=1,2,\ldots N$, 
denote a collection of dipole strengths, and $d_{j} \in \mathbb{R}^{2}$,
$j=1,2,\ldots N$, denote the corresponding dipole orientation vectors.
Let $k\in\mathbb{C}$ denote the wave number or the Helmholtz 
parameter. 

The Helmholtz FMM computes 
the potential $u(x) \in \mathbb{C}$ 
given by

.. math::
   :label: helm_nbody

    u(x) = \sum_{j=1}^{N} c_{j} H_{0}^{(1)}(k\|x-x_{j}\|) - v_{j} d_{j}\cdot \nabla H_{0}^{(1)}(k\|x-x_{j}\|)  \, , 

and its gradient $\nabla u(x) \in \mathbb{C}^{2}$
at the source and target locations, where $H_{0}^{(1)}$ is the Hankel function
of the first kind of order $0$. When $x=x_{j}$, the term
corresponding to $x_{j}$ is dropped from the sum.


Biharmonic FMM
***************
Let $c_{j} = (c_{j,1}, c_{j,2})\in \mathbb{C}^2$, 
$j=1,2,\ldots N$, 
denote a collection of charge strengths, and 
$v_{j} = (v_{j,1}, v_{j,2}, v_{j,3}) \in \mathbb{C}^{3}$,
$j=1,2,\ldots N$, 
denote a collection of dipole strengths.

The Biharmonic FMM computes 
the potential $u(x)$ and its `gradient` = 
$(P_{z} \frac{\mathrm{d}}{\mathrm{d}z}, P_{\overline{z}} \frac{\mathrm{d}}{\mathrm{d}z}, \frac{\mathrm{d}}{\mathrm{d}\overline{z}})$
given by

.. math::
   :label: biharm_nbody

    u(z) &= \sum_{j=1}^{N} c_{j,1} \log{\|z - \xi_{j}\|} + 
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
corresponding to $\xi_{j}$ is dropped from the sum. Here
$z = x(1) + i x(2)$, and $\xi_{j} = x_{j}(1) + ix_{j}(2)$, are the complex
numbers corresponding to $x$ and $x_{j}$ respectively.
When $x=x_{j}$, the term
corresponding to $x_{j}$ is dropped from the sum.
In the above, $(P_{z} \frac{\mathrm{d}}{\mathrm{d}z}$
is the part of the gradient $\frac{\mathrm{d}}{\mathrm{d} z}$
which is purely a function of $z$, and similarly
$(P_{\overline{z}} \frac{\mathrm{d}}{\mathrm{d}z}$
is the part of the gradient $\frac{\mathrm{d}}{\mathrm{d} z}$
which is purely a function of $\overline{z}$.


Modified Biharmonic FMM
************************

Let $G^{\textrm{mbh}}(x,y)$ denote the modified biharmonic
Green's function given by

.. math::
    G^{\textrm{mbh}}(x,y) = \frac{1}{2\pi \beta^2}\left(K_{0}(\beta \|x-y \|) - \log{\|x-y\|}\right)

where $K_{0}$ is the modified Bessel function of order $0$, and $\beta$ is the
Modified Biharmonic wavenumber.


Let $c_{j} \in \mathbb{R}$, 
denote a collection of charge strengths, 
$v_{j} \in \mathbb{R}$, 
denote a collection of dipole strengths,
$d_{j} = (d_{j,1}, d_{j,2})$ denote a collection
of dipole vectors, 
$q_{j} \in \mathbb{R}$ denote a collection of 
quadrupole strengths, 
$w_{j} = (w_{j,1}, w_{j,2}, w_{j,3}) \in \mathbb{R}^{3}$, 
denote a collection of quadrupole orientation vectors, 
$o_{j}$ denote a collection of octopole strengths, and
$p_{j} = (p_{j,1}, p_{j,2}, p_{j,3}, p_{j,4}) \in \mathbb{R}^{4}$, 
denote a collection of octopole strengths. For all the vectors,
$j=1,2,\ldots N$.

The Modified Biharmonic FMM computes 
the potential $u(x)\in \mathbb{R}$ 
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
corresponding to $x_{j}$ is dropped from the sum.

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

$P^{\textrm{stok}}(x,y)$ denote the associated pressure tensor

.. math::
    P(x,y) = \frac{1}{\|x-y \|^2}\begin{bmatrix}
    (x_{1}-y_{1}) &
    (x_{2}-y_{2})
    \end{bmatrix} \, ,

$T^{\textrm{stok}}(x,y)$ denote the Stresslet whose action of a vector $v$
is given by

.. math::
    v\cdot T^{\textrm{stok}}(x,y) = -\frac{2 v \cdot (x-y)}{\|x-y \|^4}
    \begin{bmatrix}
    (x_{1} - y_{1})^2 & (x_{1}-y_{1})(x_{2}-y_{2}) \\
    (x_{1}-y_{1})(x_{2}-y_{2}) & (x_{2} - y_{2})^2
    \end{bmatrix} \, ,

and finally let $\Pi^{\textrm{stok}} (x,y)$ denote its associated pressure tensor given by

.. math::
    v\cdot \Pi(x,y)^{\textrm{stok}} = -\frac{v}{\|x-y \|^2} + \frac{2 v \cdot(x-y)}{\|x-y \|^4}
    \begin{bmatrix}
    (x_{1}-y_{1}) \\
    (x_{2}-y_{2})
    \end{bmatrix}

Let $c_{j} \in \mathbb{R}^2$, 
$j=1,2,\ldots N$, 
denote a collection of Stokeslet strengths, $v_{j} \in \mathbb{R}^2$,
$j=1,2,\ldots N$, 
denote a collection of Stresslet strengths, and $d_{j} \in \mathbb{R}^{2}$,
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
corresponding to $x_{j}$ is dropped from the sum.

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
   and targets are contained in a cube with side length $L$, then
   for all $x$ such that $\| x-x_{j} \| \leq L \varepsilon_{\textrm{mach}}$,
   the term corresponding to $x_{j}$ is dropped from the sum.
   Here $\varepsilon_{\textrm{mach}} = 2^{-52}$ is machine precision.

