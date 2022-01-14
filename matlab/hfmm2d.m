function [pottarg,gradtarg] = hfmm2d(eps,zk,srcinfo,targ)
    charges = srcinfo.charges;
    dipstr = srcinfo.dipstr;
    dipvec = srcinfo.dipvec;
    sources = srcinfo.sources;
    [~,ns] = size(sources);
    [~,nt] = size(targ);
    pottarg = complex(zeros(1,nt));
    gradtarg = complex(zeros(2,nt));
    ier = 0;
    
    mex_id_ = 'hfmm2d_t_cd_g(i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[x], i dcomplex[x], i double[xx], i int[x], i double[xx], io dcomplex[x], io dcomplex[xx], io int[x])';
[pottarg, gradtarg, ier] = fmm2d(mex_id_, eps, zk, ns, sources, charges, dipstr, dipvec, nt, targ, pottarg, gradtarg, ier, 1, 1, 1, 2, ns, ns, ns, 2, ns, 1, 2, nt, nt, 2, nt, 1);
    pottarg = reshape(pottarg,[1,nt]);
    gradtarg = reshape(gradtarg, [2,nt]);

end



