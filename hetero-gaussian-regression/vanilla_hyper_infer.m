function [sigma, l, seps] = vanilla_hyper_infer(dim, obi, filledValue, sigma)
    nf = dim(1);    
    l = 50;
    seps = 0.001;
    loss = @(X)evidence(X, obi, nf, filledValue);
    X0 = log([sigma; l; seps]);
    optmethod = 1;
    checkGradient = false;
    if optmethod ==1
        options = optimoptions('fmincon','Algorithm','sqp',...
            'SpecifyObjectiveGradient',true,...
            'Display', 'off',...
            'CheckGradients', checkGradient,...
            'FiniteDifferenceType', 'central',...
            'MaxIterations', 25 );
        lb = -10*ones(size(X0)); %lower bound 
        lb(2) = 3;
        ub = 5*ones(size(X0)); %upper bound
        [x_min, fval, exitflag, output] = fmincon(loss,X0,[],[],[],[],lb, ub,[], options);
    elseif optmethod == 2
        options = optimoptions('fminunc','Algorithm','quasi-newton',...
            'SpecifyObjectiveGradient',true,...
            'Display', 'off',...
            'CheckGradients', checkGradient,...
            'FiniteDifferenceType', 'central',...
            'MaxFunctionEvaluations', 500,...
            'MaxIterations', 100 );
        [x_min, fval, exitflag, output] = fminunc(loss,X0, options);
    else
        x_min = X0;
    end
    sigma = exp(x_min(1));
    l = exp(x_min(2));
    seps = exp(x_min(3));
%     logp = -fval; 
end    

function [f, g] = evidence(X, obi, nf, filledValue)
    filled = length(filledValue);
    sigma = exp(X(1));
    l = exp(X(2));
    seps = exp(X(3));
    
    cov = zeros(filled);
    cov2 = cov;
    for i = 1:filled
        for j = 1:filled
            [x1, y1] = unvec(obi(i), nf);
            [x2, y2] = unvec(obi(j), nf);

            dist = (x1-x2)^2+(y1-y2)^2;
            cov(i,j) = sigma^2*exp(-dist/(2*l^2));
            cov2(i,j) = cov(i,j)*dist/l^3;
        end
    end
    K = cov + diag(seps^2*ones(filled,1));
    try
        L = chol(K);
    catch 
        L = chol(nearestSPD(K));
    end
    assert(norm(L'*L - K)<1e-9);
    beta = L'\filledValue;
    f = 0.5*(beta'*beta)+sum(log(diag(L))); 
    
    if ( nargout > 1 )
%         alpha = generalizedMinimalResidualMethod(K, filledValue);
        alpha = K\filledValue;
        aat = alpha*alpha';
        invK = inv(L)*inv(L)';
        ami = aat - invK;

        %calculate partial derivative for sigma
        psig = trace(ami*cov)/sigma;

        %calculate partial derivative for seps
        pseps = trace(ami)*seps;

        %partial derivative for l
        pl = 0.5*trace(ami*cov2);

        g = [-psig*sigma; -pl*l; -pseps*seps];
    end
end

