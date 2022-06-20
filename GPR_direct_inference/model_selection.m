%% initialize parameter
if nitr == 1
    sigma = 0.4;
    l = 0.0007;
    h = sqrt(1/l);
    seps = 0.1*ones(nent,1); %\sigma_\eps
    logp = -inf;
    cov = zeros(nent, nent);
    cov2 = zeros(nent, nent); %this covariance matrix is for gradient calculation
end
step_size = 0.0001;
itr = 0;
diff = [];

%% training hyperparameter
while itr<30
    %calculate K for this round of update
    for i = 1:nent
        for j = 1:nent
            [x1, y1] = unvec(i, dim(1));
            [x2, y2] = unvec(j, dim(1));
            dist = (x1-x2)^2+(y2-y1)^2;
            cov(i,j) = sigma^2*exp(-l*dist/2);
            cov2(i,j) = -sigma^2*dist*exp(-l*dist/2)/2;
        end
    end
    origcov = cov;
    cov = cov + seps^2*eye(nent);
    K = cov(obi, obi);
    L = chol(nearestSPD(K))';
    beta = L\filledValue(1:filled)';
    alpha = L'\beta;
    aat = alpha*alpha';
    invK = K\eye(filled);
    
    %calculate partial derivative for sigma
    psig = trace((aat-invK)*origcov(obi, obi))/sigma;
    
    %partial derivative for l
    pl = 0.5*trace((aat-invK)*cov2(obi,obi));
    ph = -2*pl/h^3;
    
    %partial derivative for \sigma_\eps
    pseps = seps*trace(aat-invK);
    
    %update the hyperparameters
    sigma = sigma + step_size*psig;
    h = h + step_size*ph;
    l = 1/h^2;
    seps = seps + step_size*pseps;
    
    %calculate the log probability
    oldlogp = logp;
    % do not directly calculate determinant, it will blow up to inf
    % instead use cholesky decomposition
    % should have term - nent*log(2*pi)/2, but omit since it is a constant
    logp = -0.5*(beta'*beta)-0.5*sum(log(diag(L)));
    fprintf('iter num %i, loss: %1.2e\n',itr,-logp);
    itr = itr+1;
    diff(end+1) = abs(oldlogp - logp);
    if numel(diff)>2 && mean(diff(end-2:end))< 0.01
        break
    end
    if logp == inf ||logp == -inf
        break
    end
    if itr>15
        step_size = step_size/2;
    end
end