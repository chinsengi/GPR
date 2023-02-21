load('prMatrix.mat')
phi=randperm(numel(X(:)),200);
M=zeros(size(X));
M(phi)=X(phi);

opts.MAX_ITERS=2000;
opts.GEN_PLOTS=false;
opts.QUIET = true;

l1=0.5;
g=0.06;
r1=0;
scale=1;

l2=0;
p1=size(X,1);
q1=4;
p2=2;
q2=1;
X_c = DRM(M, r1, l1,l2,g,p1,q1,p2,q2,opts);
MSE=norm(X-X_c,'fro')^2/(numel(X(:)));