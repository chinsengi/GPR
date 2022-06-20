%% initialize parameter
logp = -inf;
filled = length(obi); %number of observed entry 
cov = zeros(filled, filled);
cov2 = zeros(filled, filled); %this covariance matrix is for gradient calculation

if ~heterosig
    sigma = ones(filled,1)*sigma;
end
if ~heterol
    l = l*ones(filled,1);
end
if ~heteroseps
    seps = seps*ones(filled,1); %\sigma_\eps
end


psig = zeros(filled,1);
pseps = zeros(filled,1);
plogh = zeros(filled,1);
pl = zeros(filled, 1);
verbose = true;

step_size = 1e-5;
%% training hyperparameter
for itr = 1:600
    %calculate K for this round of update
    for i = 1:filled
        for j = 1:filled
            [x1, y1] = unvec(obi(i), nf);
            [x2, y2] = unvec(obi(j), nf);
            dist = (x1-x2)^2+(y2-y1)^2;

            cov(i,j) = sigma(i)*sigma(j)*(2*l(i)*l(j)/(l(i)^2+l(j)^2))...
                                      *exp(-dist/(l(i)^2+l(j)^2));
            cov2(i,j) = cov(i,j)*((-2*l(j)*l(i)^2 + 2*l(j)^3)/...
                 ((l(i)^2+l(j)^2)*l(i)*l(j))+4*l(i)*dist/(l(i)^2+l(j)^2)^2);
        end
    end
    trueCov = cov;
    K = cov + diag(seps.^2);
%     K = preprocess(K, false);
%     L = chol(nearestSPD(K))';
    % start generalized residual method
%     beta = L\filledValue(1:filled)';
%     alpha = L'\beta;
    alpha = generalizedMinimalResidualMethod(K, filledValue');
    aat = alpha*alpha';
%     invK = inv(L)*inv(L)';
    invK = inv(K);
    ami = aat - invK;

    %save old derivative
    oldpsig = psig;
    oldpseps = pseps;
    oldpl = pl;

    %calculate partial derivative for sigma
    if heterosig
        psig = (1./sigma).*diag(ami*trueCov);
    else 
        psig = trace(ami*trueCov)/sigma(1);
    end

    %partial derivative for \sigma_\eps
    if heteroseps
        pseps = seps.*diag(ami);
    else 
        pseps = seps*trace(ami);
    end

    %partial derivative for l
    if heterol
        for i = 1:filled
            tmp = zeros(filled, filled);
            % when i = j , the cov2 is zero anyway. 
            tmp(i,:) = cov2(i,:);
            tmp(:,i) = cov2(i,:);
            pl(i) = 0.25*trace(ami*tmp);
        end
    else
        pl = 0.5*trace(ami*cov2);
    end

    %update the hyperparameters 
    sigma = sigma + step_size*psig;
    seps = seps + step_size*pseps;
    l = l + step_size*pl;
    

    %calculate the log probability
    oldlogp = logp;
    % do not directly calculate determinant, it will blow up to inf
    % instead use cholesky decomposition
    % a constant term omitted
    L = chol(nearestSPD(K))';
    beta = L\filledValue(1:filled)';
    logp = -0.5*(beta'*beta)-sum(log(diag(L))); 
    if verbose && mod(itr, 50)==0
        fprintf('iter1 num %i, loss: %1.5e\n',itr,-logp);
    end

    %implement backtracking line search to determine stepsize
    if -logp > -oldlogp - 0.5*step_size*norm([oldpsig; oldpseps; oldpl])^2
        step_size = step_size*0.5;
    end

    if logp == inf ||logp == -inf ||step_size<1e-9
        break
    end           
end

%% learning hyperparameter
finalseps = seps;
finalsigma = sigma;
finall = l;

