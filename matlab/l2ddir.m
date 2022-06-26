function [U] = l2ddir(srcinfo,targ,pgt)
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
%  
%  The sum is evaluated directly - (slow code for testing)
% 
%  Args:
%
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
%     *  srcinfo.dipvec: double(nd,2,n) 
%           dipole orientation vectors, $v_{j}$ (optional
%           default - term corresponding to dipoles dropped) 
%  
%  -  targ: double(3,nt)
%        target locations, $t_{i}$ 
%  -  pgt: integer
%        | target eval flag 
%        | potential at targets evaluated if pgt = 1
%        | potential and gradient at targets evaluated if pgt=2  
%  
%  Returns:
%  
%  -  U.pottarg: potential at target locations, if requested, $u(t_{i})$
%  -  U.gradtarg: gradient at target locations, if requested, $\nabla u(t_{i})$
%  -  U.hesstarg: hessian at target locations, if requested, $\nabla \nabla u(t_{i})$
 

  thresh = 1e-15;
  sources = srcinfo.sources;
  [m,ns] = size(sources);
  assert(m==2,'The first dimension of sources must be 2');
  if(~isfield(srcinfo,'nd'))
    nd = 1;
  end
  if(isfield(srcinfo,'nd'))
    nd = srcinfo.nd;
  end


  pottarg = complex(zeros(nd,1));
  gradtarg = complex(zeros(nd*2,1));
  hesstarg = complex(zeros(nd*3,1));
  [m,nt] = size(targ);
  assert(m==2,'First dimension of targets must be 2');
  if(pgt >=1), pottarg = complex(zeros(nd,nt)); end;
  if(pgt >= 2), gradtarg = complex(zeros(nd*2,nt)); end;
  if(pgt == 3), hesstarg = complex(zeros(nd*3,nt)); end;

  if(pgt ==0), disp('Nothing to compute, set pgt to 1,2 or 3'); return; end;

  if(isfield(srcinfo,'charges'))
    ifcharge = 1;
    charges = srcinfo.charges;
    if(nd==1), assert(length(charges)==ns,'Charges must be same length as second dimension of sources'); end;
    if(nd>1), [a,b] = size(charges); assert(a==nd && b==ns,'Charges must be of shape [nd,ns] where nd is the number of densities, and ns is the number of sources'); end;
  else
    ifcharge = 0;
    charges = complex(zeros(nd,1));
  end

  if(isfield(srcinfo,'dipstr') || isfield(srcinfo,'dipvec'))
    ifdipole = 1;
    dipstr = srcinfo.dipstr;
    if(nd==1), assert(length(dipstr)==ns,'Dipole strength must be same length as second dimension of sources'); end;
    if(nd>1), [a,b] = size(dipstr); assert(a==nd && b==ns,'Dipstr must be of shape [nd,ns] where nd is the number of densities, and ns is the number of sources'); end;
    dipvec = srcinfo.dipvec;
    if(nd == 1), [a,b] = size(squeeze(dipvec)); assert(a==2 && b==ns,'Dipvec must be of shape[2,ns], where ns is the number of sources'); end;
    if(nd>1), [a,b,c] = size(dipvec); assert(a==nd && b==2 && c==ns, 'Dipvec must be of shape[nd,2,ns], where nd is number of densities, and ns is the number of sources'); end;
    dipvec = reshape(dipvec,[2*nd,ns]);
  else
    ifdipole = 0;
    dipvec = zeros(nd*2,ns);
    dipstr = complex(zeros(nd,ns));
  end

  nd2 = 2*nd;
  nd3 = 3*nd;

  if(pgt == 1)
    if(ifcharge==1 && ifdipole == 0)
      mex_id_ = 'l2d_directcp(i int[x], i double[xx], i int[x], i dcomplex[xx], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pottarg] = fmm2d(mex_id_, nd, sources, ns, charges, targ, nt, pottarg, thresh, 1, 2, ns, 1, nd, ns, 2, nt, 1, nd, nt, 1);
    end
    if(ifcharge==0 && ifdipole == 1)
      mex_id_ = 'l2d_directdp(i int[x], i double[xx], i int[x], i dcomplex[xx], i double[xx], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pottarg] = fmm2d(mex_id_, nd, sources, ns, dipstr, dipvec, targ, nt, pottarg, thresh, 1, 2, ns, 1, nd, ns, nd2, ns, 2, nt, 1, nd, nt, 1);
    end
    if(ifcharge==1 && ifdipole == 1)
      mex_id_ = 'l2d_directcdp(i int[x], i double[xx], i int[x], i dcomplex[xx], i dcomplex[xx], i double[xx], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pottarg] = fmm2d(mex_id_, nd, sources, ns, charges, dipstr, dipvec, targ, nt, pottarg, thresh, 1, 2, ns, 1, nd, ns, nd, ns, nd2, ns, 2, nt, 1, nd, nt, 1);
    end
    U.pottarg = pottarg;
  end
  if(pgt == 2)
    if(ifcharge==1 && ifdipole == 0)
      mex_id_ = 'l2d_directcg(i int[x], i double[xx], i int[x], i dcomplex[xx], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg] = fmm2d(mex_id_, nd, sources, ns, charges, targ, nt, pottarg, gradtarg, thresh, 1, 2, ns, 1, nd, ns, 2, nt, 1, nd, nt, nd2, nt, 1);
    end
    if(ifcharge==0 && ifdipole == 1)
      mex_id_ = 'l2d_directdg(i int[x], i double[xx], i int[x], i dcomplex[xx], i double[xx], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg] = fmm2d(mex_id_, nd, sources, ns, dipstr, dipvec, targ, nt, pottarg, gradtarg, thresh, 1, 2, ns, 1, nd, ns, nd2, ns, 2, nt, 1, nd, nt, nd2, nt, 1);
    end
    if(ifcharge==1 && ifdipole == 1)
      mex_id_ = 'l2d_directcdg(i int[x], i double[xx], i int[x], i dcomplex[xx], i dcomplex[xx], i double[xx], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg] = fmm2d(mex_id_, nd, sources, ns, charges, dipstr, dipvec, targ, nt, pottarg, gradtarg, thresh, 1, 2, ns, 1, nd, ns, nd, ns, nd2, ns, 2, nt, 1, nd, nt, nd2, nt, 1);
    end
    U.pottarg = pottarg;
    U.gradtarg = squeeze(reshape(gradtarg,[nd,2,nt]));
  end
  if(pgt == 3)
    if(ifcharge==1 && ifdipole == 0)
      mex_id_ = 'l2d_directch(i int[x], i double[xx], i int[x], i dcomplex[xx], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg, hesstarg] = fmm2d(mex_id_, nd, sources, ns, charges, targ, nt, pottarg, gradtarg, hesstarg, thresh, 1, 2, ns, 1, nd, ns, 2, nt, 1, nd, nt, nd2, nt, nd3, nt, 1);
    end
    if(ifcharge==0 && ifdipole == 1)
      mex_id_ = 'l2d_directdh(i int[x], i double[xx], i int[x], i dcomplex[xx], i double[xx], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg, hesstarg] = fmm2d(mex_id_, nd, sources, ns, dipstr, dipvec, targ, nt, pottarg, gradtarg, hesstarg, thresh, 1, 2, ns, 1, nd, ns, nd2, ns, 2, nt, 1, nd, nt, nd2, nt, nd3, nt, 1);
    end
    if(ifcharge==1 && ifdipole == 1)
      mex_id_ = 'l2d_directcdh(i int[x], i double[xx], i int[x], i dcomplex[xx], i dcomplex[xx], i double[xx], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg, hesstarg] = fmm2d(mex_id_, nd, sources, ns, charges, dipstr, dipvec, targ, nt, pottarg, gradtarg, hesstarg, thresh, 1, 2, ns, 1, nd, ns, nd, ns, nd2, ns, 2, nt, 1, nd, nt, nd2, nt, nd3, nt, 1);
    end
    U.pottarg = pottarg;
    U.gradtarg = squeeze(reshape(gradtarg,[nd,2,nt]));
    U.hesstarg = squeeze(reshape(hesstarg,[nd,3,nt]));
  end
end
%
%
% ---------------------------------------------------------------------
