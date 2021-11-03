function [pottarg,gradtarg] = h2ddir(eps,zk,srcinfo,targ)
    nd = 1;
    charges = srcinfo.charges;
    dipstr = srcinfo.dipstr;
    dipvec = srcinfo.dipvec;
    sources = srcinfo.sources;
    [~,ns] = size(sources);
    [~,nt] = size(targ);
    pottarg = complex(zeros(1,nt));
    gradtarg = complex(zeros(2,nt));
    ier = 0;
    nd2 = 2*nd;
    thresh = 1e-16;
    
    mex_id_ = 'h2d_directcdg(i int[x], i dcomplex[x], i double[xx], i int[x], i dcomplex[xx], i dcomplex[xx], i double[xx], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg] = fmm2d(mex_id_, nd, zk, sources, ns, charges, dipstr, dipvec, targ, nt, pottarg, gradtarg, thresh, 1, 1, 2, ns, 1, nd, ns, nd, ns, nd2, ns, 2, nt, 1, nd, nt, nd2, nt, 1);
    pottarg = reshape(pottarg,[1,nt]);
    gradtarg = reshape(gradtarg, [2,nt]);
    
end
