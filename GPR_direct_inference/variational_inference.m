n = 1;
nmax = 60000;
step_size = 1e-3;

executionEnvironment = "auto";
nvar = dim(1)*Rank;

while n < nmax
    [KL,gmean, gB, reconst, reg] = dlfeval(@varinf, mu, B, Rank, dim(1), dW, obi, cov, covu);
    learnRate = step_size;
    gradDecay = 0.75;
    sqGradDecay = 0.95;
    averageGrad = [];
    averageSqGrad = [];
    averageGrad1 = [];
    averageSqGrad1 = [];

    %use adam update scheme
    [mu,averageGrad,averageSqGrad] = adamupdate(mu,gmean,averageGrad,averageSqGrad,n,learnRate,gradDecay,sqGradDecay);
%     [B,averageGrad1,averageSqGrad1] = adamupdate(B,gB,averageGrad1,averageSqGrad1,n,learnRate,gradDecay,sqGradDecay);
    
    % fixed step size is unstable
%     step_size = 1e-5;
%     mu = mu - gmean*step_size;
%     B = B - gB*step_size;
    
    % project B onto Bp such that Bp*Bp' has squared exponential block
    % diagonal
%     B = extractdata(B);
%     Sigma = B*B';
%     Sigma(1:nvar, 1:nvar) = (trace(Sigma(1:nvar, 1:nvar))/nvar)*covu;
%     Sigma(nvar+1:end, nvar+1:end) = (trace(Sigma(nvar+1:end, nvar+1:end))/nvar)*covu;
%     B = dlarray(chol(nearestSPD(Sigma)));
%     

    %smooth the mean (crucial for smoothness)
    mu = extractdata(mu);
    [U, V] = convert(mu, Rank, dim(1));
    K = 0.25*ones(2);
    U(end+1,:) = U(end, :);
    U(:, end+1) = U(:, end);
    V(end+1,:) = V(end, :);
    V(:, end+1) = V(:, end);
    Usmooth = conv2(U,K,'valid');
    Vsmooth = conv2(V,K,'valid');
    mu = dlarray([reshape(Usmooth, [],1); reshape(Vsmooth, [],1)]);
    
    if mod(n, 100) == 0
        fprintf('iter num %i, KL Divergence: %1.2e, reconst: %1.2e, reg: %1.2e\n',n,KL, extractdata(reconst), extractdata(reg));
    end
    fprintf('iter num %i, KL Divergence: %1.2e, reconst: %1.2e, reg: %1.2e\n',n,KL, extractdata(reconst), extractdata(reg));
    n = n+1;    
end

%m stands for mean, sigma = BB', obi is the index of observed entries 
function [KL, gmean, gB, reconst, reg] = varinf(m, B, Rank, dim, dW, obi, cov, covu) 
    KL = 0;
    gB = 0;
    
    %calculate differential entropy Hq
%     Hq = log(abs(det(extractdata(B))));
%     KL = KL-Hq;
%     gB = inv(extractdata(B))';
    
    %monte carlo gradient estimation with reparameterization trick
    ntrial = 1;
    cov = cov(obi, obi);
    
    X0 = reshape(dW', [], 1);
    X0 = X0(obi); %observed entries
    assert(sum(X0==0)==0) %check the permutation is correct
    
    expectation = 0;
    reconst = 0;
    reg = 0;
    for trial = 1:ntrial
        % sample epsilon and replace U and V by m + B*epsilon
        epsilon = mvnrnd(zeros(size(m)), eye(size(B)), 1);
%         [U, V] = convert(m+B*epsilon', Rank, dim);
        [U, V] = convert(m, Rank, dim);
        Xflat = reshape(V*U', [], 1);
        X = Xflat(obi);
        
        % calculate  Eq[ln p (R|U,V)] = -1/2 * (X-X0)'K(X-X0)
        tmp = 0.5*(X - X0)'*inv(cov)*(X-X0);
        expectation = expectation + tmp;
        reconst = reconst+tmp;
        
        %calculate Eq[p(U)]
%         vecu = reshape(U, [],1);
%         expectation = expectation + 0.5*vecu'*pinv(covu)*vecu;
        
        %calculate Eq[p(U)]
%         vecv = reshape(V, [],1);
%         expectation = expectation + 0.5*vecv'*pinv(covu)*vecv;
        
    end
    reconst = reconst/ntrial;
    expectation = expectation/ntrial;
    reg = expectation-reconst;
    KL = KL + expectation;
    gB = gB + dlgradient(expectation, B);
    gmean = dlgradient(expectation, m);
    
%     reg = penalty;
end

function [U, V] = convert(mu, rank, dim)
    vecu = mu(1:rank*dim);
    vecv = mu(rank*dim+1:end);
    U = reshape(vecu, [], rank);
    V = reshape(vecv, [], rank);
end


function norm = calnorm(A)
    norm = sum(A.*A, 'all');
end