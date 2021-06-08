%% initialize parameter
sigma = sqrt(abs(mean(h))/(nn*cd*meanRate));
l = 25;
seps = 0.1; %\sigma_\eps
noise_mu = 1;
seps_neuron = 0.05; % seps_neuron is the noise param for the original map (instead of the affine one)

loss = @(X)grad_f(X, h, vb_index, C, rf);
X0 = log([sigma; l; seps; noise_mu; seps_neuron]);
options = optimoptions('fminunc','Algorithm','quasi-newton',...
    'SpecifyObjectiveGradient',true,...
    'Display', 'final',...
    'CheckGradients', false,...
    'FiniteDifferenceType', 'central',...
    'MaxFunctionEvaluations', 1e5);
[x_min, fval, exitflag, output] = fminunc(loss,X0,options);
[sigma, l, seps, noise_mu, seps_neuron] = extract_param(x_min);
logp = -fval;

function [f, g] = grad_f(X, h, vb_index, C, rf)
    [sigma, l, seps, noise_mu, seps_neuron] = extract_param(X);
    [K, Kl, Ksn, Ksig] = Kgen(sigma, l, seps_neuron, vb_index, C, rf);
    % calculate K for this round of update
    K = K + diag(seps^2*ones(length(K),1));
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
        
        g = [psig; pl; pseps; pnmu; psn];
        g = -g.*exp(X);
    end
end

function [sigma, l, seps, noise_mu, sn] = extract_param(X)
    [sigma, l, seps, noise_mu, sn] = feval(@(x) x{:}, num2cell(exp(X)));
end
