% model selection 
optmethod = 1;
checkGradient = false;
if model_selection_method == 1
    grad_f = @grad_f_MAP;
elseif model_selection_method == 2
    grad_f = @grad_f_cval;
else
    optmethod = 3;
end

loss = @(X)grad_f(X, h, vb_index, C, rf, heteroseps);
X0 = log([sigma; l; seps; noise_mu; seps_neuron]);
if input_noise
    loss = @(X)grad_f(X, h, vb_index, C, rf, heteroseps, postgrad);
    X0 = [X0; log(innoise)];
end

fval = -loss(X0);
if optmethod ==1
    options = optimoptions('fmincon','Algorithm','sqp',...
        'SpecifyObjectiveGradient',true,...
        'Display', 'off',...
        'CheckGradients', checkGradient,...
        'FiniteDifferenceType', 'central',...
        'MaxIterations', 25 );
    lb = -5*ones(size(X0)); %lower bound
    ub = 5*ones(size(X0)); %upper bound
    [x_min, fval, exitflag, output] = fmincon(loss,X0,[],[],[],[],lb, ub,[], options);
    [sigma, l, seps, noise_mu, seps_neuron, innoise] = extract_param(x_min, nvb, heteroseps);
elseif optmethod == 2
    options = optimoptions('fminunc','Algorithm','quasi-newton',...
        'SpecifyObjectiveGradient',true,...
        'Display', 'iter',...
        'CheckGradients', checkGradient,...
        'FiniteDifferenceType', 'central',...
        'MaxFunctionEvaluations', 500,...
        'MaxIterations', 25 );
    [x_min, fval, exitflag, output] = fminunc(loss,X0, options);
    [sigma, l, seps, noise_mu, seps_neuron, innoise] = extract_param(x_min, nvb, heteroseps);
else
    x_min = X0;
    fval = loss(x_min);
end
logp = -fval;
tmp = 1;

%gradient function for cross validation
function [f, g] = grad_f_cval(X, h, vb_index, C, rf, heteroseps, postgrad)
    nvb = length(vb_index);
    input_noise = nargin > 6;
    if input_noise
        [sigma, l, seps, noise_mu, seps_neuron, innoise] = extract_param(X, nvb, heteroseps);
        [K, Kl, Ksn, Ksig, gtg] = Kgen(sigma, l, seps_neuron, vb_index, C, rf,...
                                                         innoise, postgrad);
    else
        [sigma, l, seps, noise_mu, seps_neuron] = extract_param(X, nvb, heteroseps);
        [K, Kl, Ksn, Ksig] = Kgen(sigma, l, seps_neuron, vb_index, C, rf);
    end
    % calculate K for this round of update
    if heteroseps
        K = K + diag(seps.^2);
    else
        K = K + diag(seps^2*ones(length(K),1));
    end
    try
        L = chol(K)';
    catch
        L = chol(nearestSPD(K))';
    end
    invK = inv_chol(L);
    alpha = generalizedMinimalResidualMethod(K, h-noise_mu);
    dinvK = diag(invK);
    diff = alpha./dinvK-noise_mu;
    f = -0.5*sum(log(dinvK))+ sum(diff.^2.*dinvK)/2;  
%     f = -0.5*sum(log(dinvK))+ sum(diff.^2.*dinvK)/2;
    if ( nargout > 1 )
        %intermediate function for partial derivative calculation
        derivTmp = @(ztheta)...
            (alpha-noise_mu*dinvK).*(ztheta*alpha)...
            - 0.5*(1+(alpha.^2-(noise_mu*dinvK).^2)./dinvK).*diag(ztheta*invK);
        
        %calculate partial derivative for sigma
        zsig = invK*Ksig;
        psig = derivTmp(zsig);
        psig = sum(psig./dinvK);
        
        %partial derivative for \sigma_\eps
        if heteroseps
            error('hetero noise is not implemented for cross validation');
        else 
            zseps = 2*seps*invK;
            pseps = derivTmp(zseps);
            pseps = sum(pseps./dinvK);
        end

        %partial derivative for l
        zl = invK*Kl;
        pl = derivTmp(zl);
        pl = sum(pl./dinvK);
        
        %partial derivative for noise_mu
        pnmu = sum((alpha - noise_mu*dinvK).*(1+(sum(invK, 2))./dinvK));

        %partial derivative for seps_neuron
        zsn = invK*Ksn;
        psn = derivTmp(zsn);
        psn = sum(psn./dinvK);
        
        %partial derivative for input noise variance
        if input_noise
            error('input noise not implemented for cross validation')
            pinnoise = innoise*trace(ami*gtg);
            g = [psig; pl; pseps; pnmu; psn; pinnoise];
        else
            g = [psig; pl; pseps; pnmu; psn];
%             g = [0; 0; 0; 0; 0];
        end
        g = -g.*exp(X);
    end
end

%gradient function for MAP
function [f, g] = grad_f_MAP(X, h, vb_index, C, rf, heteroseps, postgrad)
    nvb = length(vb_index);
    input_noise = nargin > 6;
    if input_noise
        [sigma, l, seps, noise_mu, seps_neuron, innoise] = extract_param(X, nvb, heteroseps);
        [K, Kl, Ksn, Ksig, gtg] = Kgen(sigma, l, seps_neuron, vb_index, C, rf,...
                                                         innoise, postgrad);
    else
        [sigma, l, seps, noise_mu, seps_neuron] = extract_param(X, nvb, heteroseps);
        [K, Kl, Ksn, Ksig] = Kgen(sigma, l, seps_neuron, vb_index, C, rf);
    end
    % calculate K for this round of update
    if heteroseps
        K = K + diag(seps.^2);
    else
        K = K + diag(seps^2*ones(length(K),1));
    end
    if any(isnan(K(:))) || any(isinf(K(:)))
        K(isnan(K)) = 0;
        K(isinf(K)) = 0;
    end
    try
        L = chol(K)';
    catch
        L = chol(nearestSPD(K))';
    end
    beta = L\(h - noise_mu);
    theta = 1;
    f = 0.5*(beta'*beta)+sum(log(diag(L)))/theta;
    
    if ( nargout > 1 )
        invK = inv_chol(L);
        alpha = generalizedMinimalResidualMethod(K, h-noise_mu);
        aat = alpha*alpha';
        ami = aat - invK/theta;

        %calculate partial derivative for sigma
        psig = 0.5*trace(ami*Ksig);

        %partial derivative for \sigma_\eps
        if heteroseps
            pseps = seps.*diag(ami);
        else 
            pseps = seps*trace(ami);
        end

        %partial derivative for l
        pl = 0.5*trace(ami*Kl);

        %partial derivative for noise_mu
        pnmu = sum((h'-noise_mu)*invK);

        %partial derivative for seps_neuron
        psn = 0.5*trace(ami*Ksn);
        
        %partial derivative for input noise variance
        if input_noise
            pinnoise = innoise*trace(ami*gtg);
            g = [psig; pl; pseps; pnmu; psn; pinnoise];
        else
            g = [psig; pl; pseps; pnmu; psn];
%             g = [0; 0; pseps; 0; 0];
        end
        g = -g.*exp(X);
    end
end


function [sigma, l, seps, noise_mu, sn, innoise] = extract_param(X, nvb, heteroseps)
    % it's a hack, so that we don't need to know if there is input noise
    X = [X; 0];
    X = exp(X);
    if heteroseps
        sigma = X(1); l = X(2);
        seps = X(3:2+nvb);
        [noise_mu, sn, innoise] = feval(@(x) x{:}, num2cell(X(3+nvb:end)));
    else
        [sigma, l, seps, noise_mu, sn, innoise] = feval(@(x) x{:}, num2cell(X));
    end
end
