function [U] = st2ddir(srcinfo,targ,ifppregtarg)
% ST2dDIR    Direct (slow) 2d Stokes kernel sums (reference for STFMM2d).
%
%  U = st2ddir(srcinfo,targ,ifppregtarg)
%
%  Stokes direct evaluation in R^2: evaluate all pairwise particle
%  interactions with targets. This is the slow O(N^2) direct code used
%  as a reference for testing the (fast) code stfmm2d.
%
%  Kernel definitions, input and outputs arguments are identical to
%  stfmm2d (see that function for all definitions), apart from:
%  1) the first argument (eps) is absent.
%  2) there are currently no outputs at sources, meaning that U.pot, U.pre,
%     and U.grad are missing (as if ifppreg=0). In other words,
%     just targets for now, and targ is thus not an optional argument.
%
% See also: STFMM2d

  sources = srcinfo.sources;
  [m,ns] = size(sources);
  assert(m==2,'The first dimension of sources must be 2');
  if(~isfield(srcinfo,'nd'))
    nd = 1;
  end
  if(isfield(srcinfo,'nd'))
    nd = srcinfo.nd;
  end

  thresh = 1e-15;

  if( nargin <= 1 )
    return;
  else
    if( nargin <= 2 ), ifppregtarg = 3; end;
    [m,nt] = size(targ);
    assert(m==2,'First dimension of targets must be 2');
    pottarg = zeros(nd*2,nt);
    pretarg = zeros(nd,nt);
    gradtarg = zeros(nd*4,nt);
  end

  if(ifppregtarg == 0), disp('Nothing to compute, set eigher ifppregtarg to 1 or 2 or 3'); return; end;

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

  if ifstoklet == 1 && ifstrslet == 0 
    mex_id_ = 'st2ddirectstokg(i int[x], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], i double[x])';
[pottarg, pretarg, gradtarg] = fmm2d(mex_id_, nd, sources, stoklet, ns, targ, nt, pottarg, pretarg, gradtarg, thresh, 1, 2, ns, nd2, ns_stok, 1, 2, nt, 1, nd2, nt, nd, nt, nd4, nt, 1);
  else
    istress = 1;
    mex_id_ = 'st2ddirectstokstrsg(i int[x], i double[xx], i int[x], i double[xx], i int[x], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], i double[x])';
[pottarg, pretarg, gradtarg] = fmm2d(mex_id_, nd, sources, ifstoklet, stoklet, istress, strslet, strsvec, ns, targ, nt, pottarg, pretarg, gradtarg, thresh, 1, 2, ns, 1, nd2, ns_stok, 1, nd2, ns_strs, nd2, ns_strs, 1, 2, nt, 1, nd2, nt, nd, nt, nd4, nt, 1);
  end

  U.pottarg = [];
  U.pretarg = [];
  U.gradtarg = [];
  if(ifppregtarg >= 1), U.pottarg = squeeze(reshape(pottarg,[nd,2,nt])); end;
  if(ifppregtarg >= 2), U.pretarg = pretarg; end;
  if(ifppregtarg >= 3), U.gradtarg = squeeze(reshape(gradtarg,[nd,2,2,nt])); end;
end
