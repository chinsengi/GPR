%% initialize parameter
loss = @(X)evidence(X, obi, nf, filledValue, heterosig, heterol, heteroseps);  

X0 = log([sigma; l; seps]);
initial_logp = -evidence(X0, obi, nf, filledValue, heterosig, heterol, heteroseps);
optmethod = 1;
checkGradient = false;
if optmethod ==1
    options = optimoptions('fmincon','Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,...
        'Display', 'off',...
        'CheckGradients', checkGradient,...
        'FiniteDifferenceType', 'central',...
        'MaxIterations', 25 );
    lb = -5*ones(size(X0)); %lower bound 
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
    fval = -initial_logp;
end    
[sigma, l, seps] = extract_param(x_min, filled, heterosig, heterol, heteroseps);
logp = -fval;
if ~heterosig, sigma = sigma(1); end
if ~heterol, l = l(1); end
if ~heteroseps, seps = seps(1); end 
% logp = -evidence(X0, obi, nf, filledValue, heterosig, heterol, heteroseps);

%% learning hyperparameter
finalseps = seps;
finalsigma = sigma;
finall = l;

function [f, g] = evidence(X, obi, nf, filledValue, hetsig, hetl, hetseps)
    filled = length(filledValue);
    [sigma, l, seps] = extract_param(X, filled, hetsig, hetl, hetseps);
    
    cov = zeros(filled, filled);
    cov2 = zeros(filled, filled);

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
    K = cov + diag(seps.^2);
    try
        L = chol(K)';
    catch
        L = chol(nearestSPD(K))';
    end
    beta = L\filledValue';
    f = 0.5*(beta'*beta)+sum(log(diag(L))); 
    if ( nargout > 1 )
        alpha = generalizedMinimalResidualMethod(K, filledValue');
        aat = alpha*alpha';
        invK = inv(L)'*inv(L);
        ami = aat - invK;

        %calculate partial derivative for sigma
        if hetsig
            psig = (1./sigma).*diag(ami*cov);
        else 
            psig = trace(ami*cov)/sigma(1);
        end

        %partial derivative for \sigma_\eps
        if hetseps
            pseps = seps.*diag(ami);
        else 
            pseps = seps(1)*trace(ami);
        end

        %partial derivative for l
        if hetl
            pl = zeros(filled,1);
            for i = 1:filled
                tmp = zeros(filled, filled);
                % when i = j , the cov2 is zero anyway. 
                tmp(i,:) = cov2(i,:);
                tmp(:,i) = cov2(i,:)';
                % multiply by 0.25 instead of 0.5in the hetero case
                % dx^2y^2/dx = 2x^2y = 2x^3 when x = y but dx^4/dx = 4x^3
                pl(i) = 0.25*trace(ami*tmp);
            end
        else
            pl = 0.5*trace(ami*cov2);
        end
        if hetsig, psig = -psig.*sigma; else psig = -psig*sigma(1); end
        if hetl, pl = -pl.*l; else pl = -pl*l(1); end
        if hetseps, pseps = -pseps.*seps; else pseps = -pseps*seps(1); end
        
        g = [psig; pl; pseps];
    end  
end

function [sigma, l, seps] = extract_param(X, filled, hetsig, hetl, hetseps)
    cur = 1;
    if ~hetsig
        sigma = exp(X(cur))*ones(filled,1);
        cur = cur+1;
    else
        sigma = exp(X(cur:cur+filled-1));
        cur = cur+filled;
    end
    
    if ~hetl
        l = exp(X(cur))*ones(filled,1);
        cur = cur+1;
    else
        l = exp(X(cur:cur+filled-1));
        cur = cur+filled;
    end
    
    if ~hetseps
        seps = exp(X(cur))*ones(filled,1);
    else
        seps = exp(X(cur:end));
    end
end
