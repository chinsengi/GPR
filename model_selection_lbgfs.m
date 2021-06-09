loss = @(X)grad_f(X, h, vb_index, C, rf, input_noise);
X0 = log([sigma; l; seps; noise_mu; seps_neuron]);
if input_noise
    loss = @(X)grad_f(X, h, vb_index, C, rf, input_noise, postgrad);
    X0 = [X0; log(innoise)];
end

constrained = true;
if constrained
    options = optimoptions('fmincon','Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,...
        'Display', 'final',...
        'CheckGradients', false,...
        'FiniteDifferenceType', 'central');
    lb = -10*ones(size(X0));
    ub = 10*ones(size(X0));
    [x_min, fval, exitflag, output] = fmincon(loss,X0,[],[],[],[],lb, ub,[], options);
    [sigma, l, seps, noise_mu, seps_neuron, innoise] = extract_param(x_min);
else
    options = optimoptions('fminunc','Algorithm','quasi-newton',...
        'SpecifyObjectiveGradient',true,...
        'Display', 'final',...
        'CheckGradients', false,...
        'FiniteDifferenceType', 'central',...
        'MaxFunctionEvaluations', 1e5,...
        'HessUpdate', 'bfgs');
    [x_min, fval, exitflag, output] = fminunc(loss,X0, options);
    [sigma, l, seps, noise_mu, seps_neuron, innoise] = extract_param(x_min);
end
logp = -fval;

function [f, g] = grad_f(X, h, vb_index, C, rf, input_noise, postgrad)
    if input_noise
        [sigma, l, seps, noise_mu, seps_neuron, innoise] = extract_param(X);
        [K, Kl, Ksn, Ksig, gtg] = Kgen(sigma, l, seps_neuron, vb_index, C, rf,...
                                                         innoise, postgrad);
    else
        [sigma, l, seps, noise_mu, seps_neuron] = extract_param(X);
        [K, Kl, Ksn, Ksig] = Kgen(sigma, l, seps_neuron, vb_index, C, rf);
    end
    % calculate K for this round of update
    K = K + diag(seps^2*ones(length(K),1));
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
    f = 0.5*(beta'*beta)+sum(log(diag(L)));
    
    if ( nargout > 1 )
        invK = inv_chol(L);
        alpha = generalizedMinimalResidualMethod(K, h-noise_mu);
        aat = alpha*alpha';
        ami = aat - invK;

        %calculate partial derivative for sigma
        psig = 0.5*trace(ami*Ksig);

        %partial derivative for \sigma_\eps
        pseps = seps*trace(ami);

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
        end
        g = -g.*exp(X);
    end
end

function [sigma, l, seps, noise_mu, sn, innoise] = extract_param(X)
    % it's a hack, so that we don't need to know if there is input noise
    X = [X; 0];
    [sigma, l, seps, noise_mu, sn, innoise] = feval(@(x) x{:}, num2cell(exp(X)));
end
