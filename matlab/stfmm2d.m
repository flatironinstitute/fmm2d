function [U] = stfmm2d(eps,srcinfo,ifppreg,targ,ifppregtarg)
% STFMM2d    FMM in 2d for Stokes (viscous fluid hydrodynamic) kernels.
%
% U = stfmm2d(eps,srcinfo,ifppreg,targ,ifppregtarg)
%
%  Stokes FMM in R^2: evaluate all pairwise particle
%  interactions (ignoring self-interactions) and, possibly,
%  interactions with targets, to precision eps, using the fast
%  multipole method (linear scaling in number of points).
%      
%  This routine computes the sum for the velocity vector,
%
%       u_i(x) = sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
%                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k
%
%  for each requested evaluation point x and i=1,2,
%  where sigma^{(m)} is the Stokeslet force, mu^{(m)} is the
%  stresslet strength, and nu^{(m)} is the stresslet orientation
%  (note that each of these is a 2 vector per source point y^{(m)}).
%  The stokeslet kernel G and stresslet kernel T are defined below.
%  Repeated indices are taken as summed over 1,2, ie, Einstein
%  convention.
%
%  There are two options for the evaluation points $x$. These
%  can be the sources themselves (ifppreg > 0; see below) or other
%  target points of interest (ifppregtarg > 0; see below). Both options
%  can be selected in one call.
%  When $x=y_{(m)}$, the term corresponding to $y^{(m)}$ is dropped
%  from the sum.
%
%  Optionally, the associated pressure p(x) and 2x2 gradient tensor
%  grad u(x) are returned,
%
%      p(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
%          + sum_m PI_{jk}(x,y{(m)}) mu^{(m)}_j nu^{(m)}_k
%
%      grad_l u_i(x) = grad_l [sum_m G_{ij}(x,y^m) sigma^{(m)}_j
%                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k]
%
%  where the pressure stokeslet P and pressure stresslet PI are defined
%  below. Note that these two may be combined to get the stress tensor.
%
%  We use the following as the kernel definitions, noting that:
%     1) The dynamic viscosity (mu) is assumed to be 1.
%     2) All kernels are a factor 2pi larger than standard definitions
%        (this is for historical reasons).
%     Thus, in general, divide the velocity (potential or grad) outputs
%     by i2pi.mu, and pressure by 2pi, to recover standard definitions.
%
%  For a source y and target x, let r_i = x_i-y_i     (note sign)
%  and let r = sqrt(r_1^2 + r_2^2) 
%
%  The Stokeslet, G_{ij}, and its associated pressure tensor, P_j,
%  we define as
%
%     G_{ij}(x,y) = (-delta_{ij}log(r)  + r_i r_j / r^2 )/2
%     P_j(x,y) = r_j/r^2
%
%  The (Type I) stresslet, T_{ijk}, and its associated pressure
%  tensor, PI_{jk}, we define as
%   
%     T_{ijk}(x,y) = -2 r_i r_j r_k/ r^4
%     PI_{jk}(x,y) = -2 delta_{jk}/r^2 + 2 r_j r_k/r^4      
% 
%  Args:
%
%  -  eps: double   
%        precision requested
%  -  srcinfo: structure
%        structure containing the following info about sources:
%     *  srcinfo.sources: double(2,n)    
%           source locations, $x_{j}$
%     *  srcinfo.nd: integer
%           number of density vectors for the same points (optional, 
%           default - nd = 1)
%     *  srcinfo.stoklet: double(nd,2,n) 
%           Stokeslet charge strengths, $sigma_{j}$ (optional, 
%           default - term corresponding to Stokeslet charge strengths dropped)
%     *  srcinfo.strslet: double(nd,2,n) 
%           stresslet strengths, $mu_{j}$ (optional
%           default - term corresponding to stresslet strengths dropped) 
%     *  srcinfo.strsvec: double(nd,2,n) 
%           stresslet orientations, $nu_{j}$ (optional
%           default - term corresponding to stresslet orientations dropped) 
%  -  ifppreg: integer
%        | source eval flag
%        | potential at sources evaluated if ifppreg = 1
%        | potential and pressure at sources evaluated if ifppreg=2
%        | potential, pressure and gradient at sources evaluated if ifppreg=3
%  -  targ: double(2,nt)
%        target locations, $t_{i}$ (optional)
%  -  ifppregtarg: integer
%        | target eval flag (optional)
%        | potential at targets evaluated if ifppregtarg = 1
%        | potential and pressure at targets evaluated if ifppregtarg = 2 
%        | potential, pressure and gradient at targets evaluated if
%        |   ifppregtarg = 3
%  
%  Returns:
%  
%  -  U.pot: velocity at source locations if requested
%  -  U.pre: pressure at source locations if requested
%  -  U.grad: gradient of velocity at source locations if requested
%  -  U.pottarg: velocity at target locations if requested
%  -  U.pretarg: pressure at target locations if requested
%  -  U.gradtarg: gradient of velocity at target locations if requested
%
% See also: ST2dDIR

  sources = srcinfo.sources;
  [m,ns] = size(sources);
  assert(m==2,'The first dimension of sources must be 2');
  if(~isfield(srcinfo,'nd'))
    nd = 1;
  end
  if(isfield(srcinfo,'nd'))
    nd = srcinfo.nd;
  end

  pot = zeros(nd*2,1); ns_pot = 1;
  pre = zeros(nd,1); ns_pre = 1;
  grad = zeros(nd*4,1); ns_grad = 1;

  if(ifppreg >= 1), pot = zeros(nd*2,ns); ns_pot = ns; end;
  if(ifppreg >= 2), pre = zeros(nd,ns); ns_pre = ns; end;
  if(ifppreg >= 3), grad = zeros(nd*4,ns); ns_grad = ns; end;

  pottarg = zeros(nd*2,1); nt_pot = 1;
  pretarg = zeros(nd,1); nt_pre = 1;
  gradtarg = zeros(nd*4,1); nt_grad = 1;
  if( nargin <= 3 )
    nt = 0;
    ifppregtarg = 0;
    targ = zeros(2,0);
  else
    if( nargin <= 4 ), ifppregtarg = 0; end;
    [m,nt] = size(targ);
    assert(m==2,'First dimension of targets must be 2');
    if(ifppregtarg >= 1), pottarg = zeros(nd*2,nt); nt_pot = nt; end;
    if(ifppregtarg >= 2), pretarg = zeros(nd,nt); nt_pre = nt; end;
    if(ifppregtarg >= 3), gradtarg = zeros(nd*4,nt); nt_grad = nt; end;
  end

  if(ifppreg ==0 && ifppregtarg ==0), disp('Nothing to compute, set eigher ifppreg or ifppregtarg to 1 or 2 or 3'); return; end;

  if(isfield(srcinfo,'stoklet'))
    ifstoklet = 1;
    ns_stok = ns;
    stoklet = srcinfo.stoklet;
    if(nd == 1), [a,b] = size(squeeze(stoklet)); assert(a==2 && b==ns,'Stoklet must be of shape[2,ns], where ns is the number of sources'); end;
    if(nd>1), [a,b,c] = size(stoklet); assert(a==nd && b==2 && c==ns, 'Stoklet must be of shape[nd,2,ns], where nd is number of densities, and ns is the number of sources'); end;
    stoklet = reshape(stoklet,[2*nd,ns]);
  else
    ifstoklet = 0;
    ns_stok = 1;
    stoklet = zeros(nd*2,1);
  end

  if(isfield(srcinfo,'strslet') && isfield(srcinfo,'strsvec'))
    ifstrslet = 1;
    ns_strs = ns;
    strslet = srcinfo.strslet;
    strsvec = srcinfo.strsvec;
    if(nd == 1), [a,b] = size(squeeze(strslet)); assert(a==2 && b==ns,'Strslet must be of shape[2,ns], where ns is the number of sources'); end;
    if(nd == 1), [a,b] = size(squeeze(strsvec)); assert(a==2 && b==ns,'Strsvec must be of shape[2,ns], where ns is the number of sources'); end;
    if(nd>1), [a,b,c] = size(strslet); assert(a==nd && b==2 && c==ns, 'Strslet must be of shape[nd,2,ns], where nd is number of densities, and ns is the number of sources'); end;
    if(nd>1), [a,b,c] = size(strsvec); assert(a==nd && b==2 && c==ns, 'Strsvec must be of shape[nd,2,ns], where nd is number of densities, and ns is the number of sources'); end;
    strslet = reshape(strslet,[2*nd,ns]);
    strsvec = reshape(strsvec,[2*nd,ns]);
  else
    ifstrslet = 0;
    ns_strs = 1;
    strslet = zeros(nd*2,1);
    strsvec = zeros(nd*2,1);
  end

  nd2 = 2*nd;
  nd4 = 4*nd;
  ier = 0;

  mex_id_ = 'stfmm2d(i int[x], i double[x], i int[x], i double[xx], i int[x], i double[xx], i int[x], i double[xx], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], io int[x])';
[pot, pre, grad, pottarg, pretarg, gradtarg, ier] = fmm2d(mex_id_, nd, eps, ns, sources, ifstoklet, stoklet, ifstrslet, strslet, strsvec, ifppreg, pot, pre, grad, nt, targ, ifppregtarg, pottarg, pretarg, gradtarg, ier, 1, 1, 1, 2, ns, 1, nd2, ns_stok, 1, nd2, ns_strs, nd2, ns_strs, 1, nd2, ns_pot, nd, ns_pre, nd4, ns_grad, 1, 2, nt, 1, nd2, nt_pot, nd, nt_pre, nd4, nt_grad, 1);

  U.pot = [];
  U.pre = [];
  U.grad = [];
  U.pottarg = [];
  U.pretarg = [];
  U.gradtarg = [];
  if(ifppreg >= 1), U.pot = squeeze(reshape(pot,[nd,2,ns])); end;
  if(ifppreg >= 2), U.pre = pre; end;
  if(ifppreg >= 3), U.grad = squeeze(reshape(grad,[nd,2,2,ns])); end;
  if(ifppregtarg >= 1), U.pottarg = squeeze(reshape(pottarg,[nd,2,nt])); end;
  if(ifppregtarg >= 2), U.pretarg = pretarg; end;
  if(ifppregtarg >= 3), U.gradtarg = squeeze(reshape(gradtarg,[nd,2,2,nt])); end;

end

% ---------------------------------------------------------------------
