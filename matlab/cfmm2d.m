function [U,varargout] = cfmm2d(eps,srcinfo,pg,varargin)
%
%
%  This subroutine computes the N-body Laplace
%  interactions and its gradients in two dimensions where 
%  the interaction kernel is given by $\log(r)$
% 
%    u(x) = \frac{i}{4}\sum_{j=1}^{N} c_{j} \log(\|x-x_{j}\|) - d_{j} v_{j}
%    \cdot \nabla \left( \log(\|x-x_{j} \|) \right)
%
%  where $c_{j}$ are the charge densities, $d_{j}$ are the dipole
%  densities,
%  $v_{j}$ are the dipole orientation vectors, and
%  $x_{j}$ are the source locations.
%  When $\|x-x_{j}\| <= L \eps_{m}$, with $L$ being the size of the bounding
%  box of sources and targets and $\eps_{m}$ being machine precision, 
%  the term corresponding to $x_{j}$ is dropped
%  from the sum.
%
%  Note: Here gradients and hessians are complex gradients and hessians
%  as opposed to cartesian gradients and hessians, also dipvec is dropped
%
% 
%  Args:
%
%  -  eps: double   
%        precision requested
%  -  srcinfo: structure
%        structure containing sourceinfo
%     
%     *  srcinfo.sources: double(2,n)    
%           source locations, $x_{j}$
%     *  srcinfo.nd: integer
%           number of charge/dipole vectors (optional, 
%           default - nd = 1)
%     *  srcinfo.charges: complex(nd,n) 
%           charge densities, $c_{j}$ (optional, 
%           default - term corresponding to charges dropped)
%     *  srcinfo.dipstr: complex(nd,n)
%           dipole densities, $d_{j}$ (optional, 
%           default - term corresponding to dipoles dropped)
%  -  pg: integer
%        | source eval flag
%        | potential at sources evaluated if pg = 1
%        | potential and gradient at sources evaluated if pg=2
%        | potential, gradient, and hessians at sources evaluated if pg=3
%        
%  Optional args
%  -  targ: double(2,nt)
%        target locations, $t_{i}$ 
%  -  pgt: integer
%        | target eval flag 
%        | potential at targets evaluated if pgt = 1
%        | potential and gradient at targets evaluated if pgt=2 
%  -  opts: options structure, values in brackets indicate default
%           values wherever applicable
%        opts.ndiv: set number of points for subdivision criterion
%        opts.idivflag: set subdivision criterion (0)
%           opts.idivflag = 0, subdivide on sources only
%           opts.idivflag = 1, subdivide on targets only
%           opts.idivflag = 2, subdivide on sources and targets
%        opts.ifnear: include near (list 1) interactions (true)
%
%  
%  Returns:
%  
%  -  U.pot: potential at source locations, if requested, $u(x_{j})$
%  -  U.grad: gradient at source locations, if requested, $d/dz u(x_{j})$
%  -  U.hess: hessian at source locations, if requested, $d/dz d/dz u(x_{j})$
%  -  U.pottarg: potential at target locations, if requested, $u(t_{i})$
%  -  U.gradtarg: gradient at target locations, if requested, $d/dz u(t_{i})$
%  -  U.hesstarg: hessian at target locations, if requested, $d/dz d/dz u(t_{i})$
%
%  - ier: error code for fmm run
%  - timeinfo: time taken in each step of the fmm
%       timeinfo(1): form multipole step
%       timeinfo(2): multipole->multipole translation step
%       timeinfo(3): multipole to local translation, form local + multipole eval step
%       timeinfo(4): local->local translation step
%       timeinfo(5): local eval step
%       timeinfo(6): direct evaluation step
%
%
%  Examples:
%  U = cfmm2d(eps,srcinfo,pg)
%     Call the FMM for sources only with default arguments
%  U = cfmm2d(eps,srcinfo,pg,targ,pgt)
%     Call the FMM for sources + targets with default arguments
%  U = cfmm2d(eps,srcinfo,pg,opts)
%     Call the FMM for sources only with user specified arguments
%  U = cfmm2d(eps,srcinfo,pg,targ,pgt)
%     Call the FMM for sources + targets with user specified arguments 
%  [U,ier] = cfmm2d(eps,srcinfo,pg)
%     Call the FMM for sources only with default arguments and returns
%     the error code for the FMM as well
%  [U,ier,timeinfo] = cfmm2d(eps,srcinfo,pg)
%     Call the FMM for sources only with default arguments, returns
%     the error code for the FMM as well and the time split
%      
 


  sources = srcinfo.sources;
  [m,ns] = size(sources);
  assert(m==2,'The first dimension of sources must be 2');
  if(~isfield(srcinfo,'nd'))
    nd = 1;
  end
  if(isfield(srcinfo,'nd'))
    nd = srcinfo.nd;
  end

  pot = complex(zeros(nd,ns)); 
  grad = complex(zeros(nd,ns));
  hess = complex(zeros(nd,ns));
  
  if( nargin < 3)
    disp('Not enough input arguments, exiting\n');
    return;
  end
  if( nargin == 3 )
    nt = 0;
    pgt = 0;
    targ = zeros(2,1);
    opts = [];
  elseif (nargin == 4)
    nt = 0;
    pgt = 0;
    targ = zeros(2,1);
    opts = varargin{1};
  elseif (nargin == 5)
    targ = varargin{1};
    pgt = varargin{2};
    [m,nt] = size(targ);
    assert(m==2,'First dimension of targets must be 2');
    opts = [];
  elseif (nargin == 6)
    targ = varargin{1};
    pgt = varargin{2};
    [m,nt] = size(targ);
    assert(m==2,'First dimension of targets must be 2');
    opts = varargin{3};
  end
  ntuse = max(nt,1);
  pottarg = complex(zeros(nd,ntuse));
  gradtarg = complex(zeros(nd,ntuse));
  hesstarg = complex(zeros(nd,ntuse));


  if((pg ==0 && pgt ==0) || (ns == 0)), disp('Nothing to compute, set eigher pg or pgt to 1 or 2'); return; end;

  if(isfield(srcinfo,'charges'))
    ifcharge = 1;
    charges = srcinfo.charges;
    if(nd==1), assert(length(charges)==ns,'Charges must be same length as second dimension of sources'); end;
    if(nd>1), [a,b] = size(charges); assert(a==nd && b==ns,'Charges must be of shape [nd,ns] where nd is the number of densities, and ns is the number of sources'); end;
  else
    ifcharge = 0;
    charges = complex(zeros(nd,ns));
  end

  if(isfield(srcinfo,'dipstr'))
    ifdipole = 1;
    dipstr = srcinfo.dipstr;
    if(nd==1), assert(length(dipstr)==ns,'Dipole strength must be same length as second dimension of sources'); end;
    if(nd>1), [a,b] = size(dipstr); assert(a==nd && b==ns,'Dipstr must be of shape [nd,ns] where nd is the number of densities, and ns is the number of sources'); end;
  else
    ifdipole = 0;
    dipstr = complex(zeros(nd,ns));
  end

  ier = 0;

  ndiv = 20;
  idivflag = 0;
  mex_id_ = 'hndiv2d(i double[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], io int[x], io int[x])';
[ndiv, idivflag] = fmm2d(mex_id_, eps, ns, nt, ifcharge, ifdipole, pg, pgt, ndiv, idivflag, 1, 1, 1, 1, 1, 1, 1, 1, 1);
  if(isfield(opts,'ndiv'))
    ndiv = opts.ndiv;
  end

  if(isfield(opts,'idivflag'))
    idivflag = opts.idivflag;
  end

  ifnear = 1;
  if(isfield(opts,'ifnear'))
    ifnear = opts.ifnear;
  end
  iper = 1;
  timeinfo = zeros(8,1);
  mex_id_ = 'cfmm2d_ndiv(i int[x], i double[x], i int[x], i double[xx], i int[x], i dcomplex[xx], i int[x], i dcomplex[xx], i int[x], i int[x], io dcomplex[xx], io dcomplex[xx], io dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], io dcomplex[xx], i int[x], i int[x], i int[x], io double[x], io int[x])';
[pot, grad, hess, pottarg, gradtarg, hesstarg, timeinfo, ier] = fmm2d(mex_id_, nd, eps, ns, sources, ifcharge, charges, ifdipole, dipstr, iper, pg, pot, grad, hess, nt, targ, pgt, pottarg, gradtarg, hesstarg, ndiv, idivflag, ifnear, timeinfo, ier, 1, 1, 1, 2, ns, 1, nd, ns, 1, nd, ns, 1, 1, nd, ns, nd, ns, nd, ns, 1, 2, ntuse, 1, nd, ntuse, nd, ntuse, nd, ntuse, 1, 1, 1, 8, 1);

  U.pot = [];
  U.grad = [];
  U.hess = [];
  U.pottarg = [];
  U.gradtarg = [];
  U.hesstarg = [];
  if(pg >= 1), U.pot = squeeze(reshape(real(pot),[nd,ns])); end;
  if(pg >= 2), U.grad = squeeze(reshape(grad,[nd,ns])); end;
  if(pg >= 3), U.hess = squeeze(reshape(hess,[nd,ns])); end;
  if(pgt >= 1), U.pottarg = squeeze(reshape(real(pottarg),[nd,nt])); end;
  if(pgt >= 2), U.gradtarg = squeeze(reshape(gradtarg,[nd,nt])); end;
  if(pgt >= 3), U.hesstarg = squeeze(reshape(hesstarg,[nd,nt])); end;

  varargout{1} = ier;
  varargout{2} = timeinfo;
end

% ---------------------------------------------------------------------
