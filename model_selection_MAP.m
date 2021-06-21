loss = @(X)grad_f(X, h, vb_index, C, rf, heteroseps);
X0 = log([sigma; l; seps; noise_mu; seps_neuron]);
if input_noise
    loss = @(X)grad_f(X, h, vb_index, C, rf, heteroseps, postgrad);
    X0 = [X0; log(innoise)];
end

fval = -grad_f(X0, h, vb_index, C, rf, heteroseps)
optmethod = 3;
if optmethod ==1
    options = optimoptions('fmincon','Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,...
        'Display', 'iter',...
        'CheckGradients', false,...
        'FiniteDifferenceType', 'central');
    lb = -10*ones(size(X0));
    ub = 10*ones(size(X0));
    [x_min, fval, exitflag, output] = fmincon(loss,X0,[],[],[],[],lb, ub,[], options);
    [sigma, l, seps, noise_mu, seps_neuron, innoise] = extract_param(x_min, nvb, heteroseps);
elseif optmethod == 2
    options = optimoptions('fminunc','Algorithm','quasi-newton',...
        'SpecifyObjectiveGradient',true,...
        'Display', 'off',...
        'CheckGradients', false,...
        'FiniteDifferenceType', 'central',...
        'MaxFunctionEvaluations', 1e5,...
        'HessUpdate', 'bfgs');
    try
        [x_min, fval, exitflag, output] = fminunc(loss,X0, options);
    catch
        x_min = X0;
    end
    [sigma, l, seps, noise_mu, seps_neuron, innoise] = extract_param(x_min, nvb, heteroseps);
else
    fval = grad_f(X0, h, vb_index, C, rf, heteroseps);
end
logp = -fval
tmp = 1;

function [f, g] = grad_f(X, h, vb_index, C, rf, heteroseps, postgrad)
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
