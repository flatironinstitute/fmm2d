function [pottarg] = hfmm2d(eps,zk,srcinfo,targ)
    charges = srcinfo.charges;
    dipstr = srcinfo.dipstr;
    dipvec = srcinfo.dipvec;
    sources = srcinfo.sources;
    [~,ns] = size(sources);
    [~,nt] = size(targ);
    pottarg = complex(zeros(nt,1));
    ier = 0;
    mex_id_ = 'hfmm2d_t_cd_p(i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[x], i dcomplex[x], i double[xx], i int[x], i double[xx], io dcomplex[x], io int[x])';
[pottarg, ier] = fmm2d(mex_id_, eps, zk, ns, sources, charges, dipstr, dipvec, nt, targ, pottarg, ier, 1, 1, 1, 2, ns, ns, ns, 2, ns, 1, 2, nt, nt, 1);

end
