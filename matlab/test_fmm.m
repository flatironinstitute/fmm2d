clear srcinfo

ntests = 2;
ipass = zeros(ntests,1);
errs = zeros(ntests,1);

ns = 4000;
nt = 3999;
srcinfo.sources = rand(2,ns);
srcinfo.charges = rand(1,ns) + 1i*rand(1,ns);
srcinfo.dipstr = rand(1,ns) + 1i*rand(1,ns);
srcinfo.dipvec = rand(2,ns);

targ = rand(2,nt);
eps = 1e-5;
zk = complex(1.1);
[pottarg,gradtarg] = hfmm2d(eps,zk,srcinfo,targ);

ntest = 10;
ttmp = targ(:,1:ntest);

[pottest,gradtest] = h2ddir(eps,zk,srcinfo,ttmp);
erra = norm(pottarg(1:ntest)-pottest)^2 + norm(gradtarg(:,1:ntest)-gradtest)^2;
ra = norm(pottest)^2  +norm(gradtest)^2;
errs(1) = sqrt(erra/ra);

assert(errs(1)<eps,'Failed helmholtz test');
ipass(1) = 1;

[pottarg] = cfmm2d(eps,srcinfo,targ);
pottest = c2ddir(eps,srcinfo,ttmp);
errs(2) = norm(pottarg(1:ntest)-pottest)/norm(pottest);

assert(errs(2)<eps,'Failed helmholtz test');
ipass(2) = 1;


isum = sum(ipass);
fprintf("Successfully cleared %d out of 2 tests in fmm2d testing suite\n",isum);


