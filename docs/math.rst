Definitions 
===========
Let $x_{j} \in \mathbb{R}^{2}$, $i=1,2,\ldots N$, denote a collection
of source locations and let $t_{i} \in \mathbb{R}^{2}$ denote a collection
of target locations. 


Laplace FMM
*************
The Laplace FMM comes in three varieties

* rfmm2d: Charges, dipole strengths, potentials, their
  gradients, and hessians are all double precision 
* lfmm2d: All of the above quantities are double complex
* cfmm2d: Same as lfmm2d except that there are no dipole orientation vectors,
  and the dipole kernel is given by $1/(zt_{i} - zx_{j})$ where $zt_{i}, zx_{j}$
  are the complex versions of the source and target locations given by
  $zt_{i} = t_{i}(1) + i*t_{i}(2)$, and $zx_{j} = x_{j}(1) + i*y_{j}(2)$. 
  Moreover, the gradients, and hessians
  in this case are derivatives with respect to $z$, and in a slight abuse of
  notation, we set $d/dz log(|z|) = 1/z$
 

Laplace FMM (rfmm2d)
---------------------
Let $c_{j} \in \mathbb{R}$, 
$j=1,2,\ldots N$, 
denote a collection of charge strengths, $v_{j} \in \mathbb{R}$,
$j=1,2,\ldots N$, 
denote a collection of dipole strengths, and $d_{j} \in \mathbb{R}^{2}$,
$j=1,2,\ldots N$, denote the corresponding dipole orientation vectors.

The Laplace FMM (rfmm2d) computes 
the potential $u(x) \in \mathbb{R}$ and the its gradient $\nabla u(x) \in
\mathbb{R}^{2}$
given by

.. math::
   :label: rlap_nbody

    u(x) = \sum_{j=1}^{N} c_{j} log(\|x-x_{j}\|) - v_{j} d_{j}\cdot \nabla log(\|x-x_{j}\|)  \, , 

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
the potential $u(x) \in \mathbb{C}$ and the its gradient $\nabla u(x) \in
\mathbb{C}^{2}$
given by

.. math::
   :label: llap_nbody

    u(x) = \sum_{j=1}^{N} c_{j} log(\|x-x_{j}\|) - v_{j} d_{j}\cdot \nabla log(\|x-x_{j}\|)  \, , 

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
the potential $u(z) \in \mathbb{C}$ and the its gradient $d/dz u(z) \in
\mathbb{C}$
given by

.. math::
   :label: clap_nbody

    u(z) = \sum_{j=1}^{N} c_{j} log(\|z-zx_{j}\|) - \frac{v_{j}}{z-zx_{j}}  \, , 

at the source and target locations. When $z=zx_{j}$, the term
corresponding to $zx_{j}$ is dropped from the sum. Here
$z = x(1) + i x(2)$, and $zx_{j} = x_{j}(1) + ix_{j}(2)$, are the complex
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
the potential $u(x)$ and the its gradient $\nabla u(x)$
given by

.. math::
   :label: helm_nbody

    u(x) = \sum_{j=1}^{N} c_{j} H_{0}^{(1)}(k\|x-x_{j}\|) - v_{j} d_{j}\cdot \nabla H_{0}^{(1)}(k\|x-x_{j}\|)  \, , 

at the source and target locations, where $H_{0}^{(1)}$ is the Hankel function
of the first kind of order $0$. When $x=x_{j}$, the term
corresponding to $x_{j}$ is dropped from the sum.


Biharmonic FMM
***************

Modified Biharmonic FMM
************************

Stokes FMM
************

Vectorized versions   
*******************
The vectorized versions of the Helmholtz FMM, 
computes repeated FMMs for new charge and dipole strengths
located at the same source locations, where the potential and its
gradient are evaluated at the same set of target locations.

For example, let $c_{\ell,j}\in\mathbb{C}$, 
$j=1,2,\ldots N$, $\ell=1,2,\ldots n_{d}$
denote a collection of $n_{d}$ charge strengths, and
let $v_{\ell,j} \in \mathbb{C}$, $d_{\ell,j} \in \mathbb{R}^2$ 
denote a collection of $n_{d}$ dipole strengths and orientation vectors. 
Then the vectorized Helmholtz FMM computes the potentials $u_{\ell}(x)$ 
and its gradients $\nabla u_{\ell}(x)$ defined by the formula

.. math::
    :label: helm_nbody_vec

    u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} H_{0}^{(1)}(k\|x-x_{j}\|) - v_{\ell,j} d_{j}\cdot \nabla H_{0}^{(1)}(k\|x-x_{j}\|)  \, , 

at the source and target locations. 

.. note::

   In double precision arithmetic, two numbers which are
   within machine precision of each other cannot be
   distinguished. In order to account for this, suppose that the sources
   and targets are contained in a cube with side length $L$, then
   for all $x$ such that $\| x-x_{j} \| \leq L \varepsilon_{\textrm{mach}}$,
   the term corresponding to $x_{j}$ is dropped from the sum.
   Here $\varepsilon_{\textrm{mach}} = 2^{-52}$ is machine precision.

