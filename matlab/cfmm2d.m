function [pottarg] = cfmm2d(eps,srcinfo,targ)
    dipstr = srcinfo.dipstr;
    sources = srcinfo.sources;
    [~,ns] = size(sources);
    [~,nt] = size(targ);
    ier = 0;
    pottarg = complex(zeros(1,nt));

    mex_id_ = 'cfmm2d_t_d_p(i double[x], i int[x], i double[xx], i dcomplex[x], i int[x], i double[xx], io dcomplex[x], io int[x])';
[pottarg, ier] = fmm2d(mex_id_, eps, ns, sources, dipstr, nt, targ, pottarg, ier, 1, 1, 2, ns, ns, 1, 2, nt, nt, 1);
    pottarg = reshape(pottarg,[1,nt]);
end



