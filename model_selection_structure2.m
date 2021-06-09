%% initialize parameter
sigma = 1;
l = 0.0007;
hprime = sqrt(1/l);
seps = 0.2*ones(nent,1); %\sigma_\eps
logp = -inf;
cov = zeros(nent, nent);
cov2 = zeros(nent, nent); %this matrix is for gradient calculation

diff = [];
psig = 0;
pseps = 0;
ph = 0;

step_size = 0.005;
itr = 0;
%% training hyperparameter
while itr<100
    % calculate K for this round of update
    % magic mex
    [cov, cov2] = covgen_mex(sigma, l, cov, cov2);
    truecov = cov;
    cov = cov + diag(seps.^2);
    K = A*cov*A';
    trueK = A*truecov*A';
    K2 = A*cov2*A';
%     K = preprocess(K);
    L = chol(nearestSPD(K))';
    beta = L\h;
    alpha = L'\beta;
    aat = alpha*alpha';
    invL = inv(L);
    invK = invL'*invL;
    ami = aat - invK;

    %save old derivative
    oldpsig = psig;
    oldpseps = pseps;
    oldph = ph;

    %calculate partial derivative for sigma
    psig = trace(ami*trueK)/sigma;

    %partial derivative for \sigma_\eps
    pseps = zeros(nent, 1);
    for i = 1:nent
        tmp = A(:,i)*A(:,i)';
        pseps(i) = seps(i)*trace(ami*tmp);
    end

    %partial derivative for l
    pl = 0.5*trace(ami*K2);
    ph = -2*pl/hprime^3;

    %update the hyperparameters 
    sigma = sigma + step_size*psig;
%     seps = seps + step_size*pseps;
    hprime = hprime + step_size*ph;
    l = 1/hprime^2;

    %calculate the log probability
    oldlogp = logp;
    % do not directly calculate determinant, it will blow up to inf
    % instead use cholesky decomposition
    % should have term - nent*log(2*pi)/2, but omit since it is a constant
    logp = -0.5*(beta'*beta)-0.5*sum(log(diag(L)));    
    fprintf('iter1 num %i, loss: %1.5e\n',itr,-logp);

    %implement backtracking line search to determine stepsize
    if -logp > -oldlogp - 0.5*step_size*norm([oldpsig; oldpseps; oldph])^2
        step_size = step_size*0.5;
    end

    itr = itr+1;
    if logp == inf ||logp == -inf || step_size<1e-9
        break
    end
end

%% learning variance
%fill in learned variance
variance = zeros(dim);
for i = 1:nent
    [x,y] = unvec(i, nf);
    variance(x,y) = seps(i);
end

%% update covariance matrix
% note that we enumerate the entries by row, not column, caution when
% reshape
cov = cov+diag(reshape(variance', nent,1).^2);
% cov = cov+diag(reshape(dWvar(1:2:100, 1:2:100), nent,1));
