function L = DRM(M, rho, lambda1,lambda2,gamma,p1,q1,p2,q2,opts)
%   this reconstructs a low rank matrix   
%   given noisy incomplete measurements, solves: 
%
%   min_X (1/2) sum_ij (L_ij - M_ij)^2 + rho * norm_nuc(L)
%
%   where the sum is over the known entries in M 
%   (the non-zero entries in this simple example)
%   Parameters:
%       M - the matrix with missing value replaced with 0
%       h - the input of neurons
%       C - random connection matrix
    opts.M = M;
    opts.rho = rho;
    opts.gamma=gamma;
    opts.lambda1=lambda1;
    opts.lambda2=lambda2;
    opts.p1=p1;
    opts.q1=q1;
    opts.p2=p2;
    opts.q2=q2;
    tic
    l = apg(@grad_f, @svd_shrink, size(M,1)*size(M,2), opts);
    L = reshape(l,size(M));
end

function g = grad_f(x, opts)
    L = reshape(x,size(opts.M));
    G = (L-opts.M).*(opts.M ~= 0);
    g = opts.gamma*G(:)+opts.lambda1*blockgradient(opts.p1,opts.q1,L,opts.lambda1)+opts.lambda2*blockgradient(opts.p2,opts.q2,L,opts.lambda2);
end

function v = svd_shrink(x, t, opts)
    L = reshape(x,size(opts.M));
    [U,S,V] = svd(L,'econ');
    s = diag(S);
    S = diag(sign(s) .* max(abs(s) - t*opts.rho,0));
    L = U*S*V';
    v = L(:);
end