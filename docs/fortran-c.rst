.. _fcexmp:

Fortran interfaces
========================

-  :ref:`helm`


.. _helm:

Helmholtz FMM
--------------

The Helmholtz FMM evaluates the following potential, its
gradient and its hessian

.. math::
 

    u(x) = \sum_{j=1}^{N} c_{j} H_{0}^{(1)}(\|x-x_{j}\|) - v_{j} d_{j}\cdot \nabla H_{0}^{(1)}(\|x-x_{j}\|)  \, . 
     
Here $x_{j}$ are the source locations, $c_{j}$ are the 
charge strengths, $v_{j}$ are the dipole strengths, and
$d_{j}$ are the dipole orientation vectors.
The collection of $x$ at which the potential
and its gradient are evaluated are referred to as the
evalution points.
     

There are 27 different Fortran wrappers for the Helmholtz FMM 
to account for collection of evaluation points (sources only, 
targets only, sources+targets), interaction kernel (charges only,
dipoles only, charges + dipoles), output request (potential,
potential+gradient,potential+gradient+hessian).

For example, the subroutine to evaluate the potential and gradient, at a collection
of targets $t_{i}$ due to a collection of charges is::

   hfmm2d_t_c_g

In general, the subroutine names take the following form::

   hfmm2d_<eval-pts>_<int-ker>_<out>



- <eval-pts>: evaluation points. Collection of `x` where $u$ and its gradient is to be evaluated

    - s: Evaluate $u$ and its gradient at the source locations $x_{i}$ 
    - t: Evaluate $u$ and its gradient at $t_{i}$, a collection of target locations specified  by the user.
    - st: Evaluate $u$ and its gradient at both source and target locations $x_{i}$ and $t_{i}$.

- <int-ker>: kernel of interaction (type of sources present)

    - c: charges
    - d: dipoles
    - cd: charges + dipoles
 
- <out>: flag for evaluating potential or potential + gradient

    - p: on output only $u$ is evaluated
    - g: on output both $u$ and its gradient are evaluated
    - h: on output $u$, its gradient and its hessian are evaluated

These are all the single density routines. To get a vectorized version 
of any of the routines use::

   <subroutine name>_vec

.. note::

   For the vectorized subroutines, the charge strengths, dipole
   strengths, potentials, and gradients are interleaved as opposed to
   provided in a sequential manner. For example for three sets of charge
   strengths, they should be stored as $c_{1,1}, c_{2,1}, c_{3,1},
   c_{1,2}, c_{2,2},c_{3,2} \ldots c_{1,N}, c_{2,N}, c_{3,N}$. 


Example drivers:

-   ``examples/hfmm2d_example.f``. The corresponding makefile is
    ``examples/hfmm2d_example.make``
-   ``examples/hfmm2d_vec_example.f``. The corresponding makefile is
    ``examples/hfmm2d_vec_example.make``

   
.. container:: rttext

  `Back to top <fortran-c.html#fcexmp>`__


List of interfaces
******************
 
- Evaluation points: Sources  
                                        
  - Interaction Type: Charges           
                                        
    - Potential (:ref:`hscp`)          
    - Gradient  (:ref:`hscg`)          
    - Hessian   (:ref:`hsch`)          
                                       
  - Interaction Type: Dipoles           
                                       
    - Potential (:ref:`hsdp`)          
    - Gradient  (:ref:`hsdg`)          
    - Hessian   (:ref:`hsdh`)          
                                       
  - Interaction Type: Charges + Dipoles 
                                       
    - Potential (:ref:`hscdp`)         
    - Gradient  (:ref:`hscdg`)         
    - Hessian   (:ref:`hscdh`)          


- Evaluation points: Targets  
                                        
  - Interaction Type: Charges           
                                        
    - Potential (:ref:`htcp`)          
    - Gradient  (:ref:`htcg`)          
    - Hessian   (:ref:`htch`)          
                                       
  - Interaction Type: Dipoles           
                                       
    - Potential (:ref:`htdp`)          
    - Gradient  (:ref:`htdg`)          
    - Hessian   (:ref:`htdh`)          
                                       
  - Interaction Type: Charges + Dipoles 
                                       
    - Potential (:ref:`htcdp`)         
    - Gradient  (:ref:`htcdg`)         
    - Hessian   (:ref:`htcdh`)          

- Evaluation points: Sources + Targets  
                                        
  - Interaction Type: Charges           
                                        
    - Potential (:ref:`hstcp`)          
    - Gradient  (:ref:`hstcg`)          
    - Hessian   (:ref:`hstch`)          
                                       
  - Interaction Type: Dipoles           
                                       
    - Potential (:ref:`hstdp`)          
    - Gradient  (:ref:`hstdg`)          
    - Hessian   (:ref:`hstdh`)          
                                       
  - Interaction Type: Charges + Dipoles 
                                       
    - Potential (:ref:`hstcdp`)         
    - Gradient  (:ref:`hstcdg`)         
    - Hessian   (:ref:`hstcdh`)          


.. container:: rttext

  `Back to top <fortran-c.html#fcexmp>`__


.. include:: fortrandocs_helm.raw


