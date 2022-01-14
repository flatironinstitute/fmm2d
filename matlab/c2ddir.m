function [pottarg] = c2ddir(eps,srcinfo,targ)
    nd = 1;
    dipstr = srcinfo.dipstr;
    sources = srcinfo.sources;
    [~,ns] = size(sources);
    [~,nt] = size(targ);
    ier = 0;
    thresh = 1e-16;
    pottarg = complex(zeros(1,nt));
    mex_id_ = 'c2d_directdp(i int[x], i double[xx], i int[x], i dcomplex[xx], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pottarg] = fmm2d(mex_id_, nd, sources, ns, dipstr, targ, nt, pottarg, thresh, 1, 2, ns, 1, nd, ns, 2, nt, 1, nd, nt, 1);
    pottarg = reshape(pottarg,[1,nt]);
     
end



