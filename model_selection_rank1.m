%% initialize parameter
sigma = sqrt(abs(mean(h))/(nn*cd*meanRate));
l = 7e-4;
hprime = sqrt(1/l); %name h is used for the data
seps = 0.1*ones(nvb,1); %\sigma_\eps
noise_mu = 1;
seps_neuron = 0*ones(nf,nf);
logp = -inf;

ALPHA = 1.1;
BETA = 0.5;
USE_GRA = true;
USE_RESTART = true;
FIXED_STEP_SIZE = false;
x = cell(4,1);
x(1) = {sigma}; x(2) = {seps}; x(3) = {hprime}; x(4) = {noise_mu}; x(5) = {seps_neuron};
% step_size initialization
% if ~FIXED_STEP_SIZE
%     % Barzilai-Borwein step-size initialization:
%     [K, K2] = Kgen(sigma, l, vb_index, C, rf);
%     [psig, pseps, phprime, pnmu, invK, L] = grad_f(sigma, hprime, seps, noise_mu, h, K, K2, heteroseps);
%     g = [psig; pseps; phprime; pnmu];
%     logp = -0.5*((h-noise_mu)'*invK*(h-noise_mu))-sum(log(diag(L))); 
%     step_size = 1 / norm(g);
%     x_hat = cell(4,1);
%     x_hat(1) = {sigma + step_size*psig}; x_hat(2) = {seps + step_size*pseps};
%     x_hat(3) = {hprime + step_size*phprime}; x_hat(4) = {noise_mu + step_size*pnmu};
%     sigma = x_hat{1}; seps = x_hat{2}; hprime = x_hat{3}; noise_mu = x_hat{4};
%     x_hat = [x_hat{1}; x_hat{2}; x_hat{3}; x_hat{4}];
%     [psig, pseps, phprime, pnmu] = grad_f(sigma, hprime, seps, noise_mu, h, K, K2, heteroseps);
%     g_hat = [psig; pseps; phprime; pnmu];
%     xv = [x{1}; x{2}; x{3}; x{4}];
%     step_size = abs(( xv - x_hat )'*(g - g_hat) / norm(g - g_hat)^2);
% else 
%     step_size = 1e-6;
% end
%% training hyperparameter
for cycle = 1:1
    psig = 0; pseps = 0; phprime = 0; pnmu = 0; psn = 0;
    step_size = 0.00001;
    brake = false;
    for itr = 1:300  
        old_x = x;
        sigma = x{1}; seps = x{2}; hprime = x{3}; 
        noise_mu = x{4}; seps_neuron = x{5};
        % update the hyperparameters 
        sigma = sigma + step_size*psig;
        seps = seps + step_size*pseps;
        hprime = hprime + step_size*phprime;
        l = 1/hprime^2;
        noise_mu = noise_mu + step_size*pnmu;
%         seps_neuron = seps_neuron + step_size*psn;
        
        x(1) = {sigma}; x(2) = {seps}; x(3) = {hprime};
        x(4) = {noise_mu}; x(5) = {seps_neuron};
        
        % gradient calculation
        [K, Kl, Ksn] = Kgen(sigma, l, vb_index, C, rf, seps_neuron, nf);
        [psig, pseps, phprime, pnmu, psn, invK, L, K] = grad_f(sigma, hprime, seps, noise_mu, h, K, Kl, Ksn, heteroseps);

        %calculate the log probability
        if ~brake
            oldlogp = logp;
        end
        % do not directly calculate determinant, it will blow up to inf
        % instead use cholesky decomposition
        % should have term - nent*log(2*pi)/2, but omit since it is a constant
        logp = -0.5*((h-noise_mu)'*invK*(h-noise_mu))-sum(log(diag(L)));  
        fprintf('iter1 num %i, loss: %1.5e\n',itr,-logp);

        %save old derivative
        oldpsig = psig;
        oldpseps = pseps;
        oldph = phprime;
        oldpnmu = pnmu;
        oldpsn = psn(:);

        if (~FIXED_STEP_SIZE)      
            %implement backtracking line search to determine stepsize
            if -logp > -oldlogp% - 0.5*step_size*norm([oldpsig; oldpseps; oldph; oldpnmu])^2
                step_size = step_size*BETA;
%                 brake = true;
%                 x = old_x;
%                 if -oldlogp < 0
%                     step_size = 1e-6;
%                 end
            else
                brake = false;
%                 step_size = step_size*ALPHA;
            end
        end


        if logp == inf ||logp == -inf || step_size<1e-9
            break;
        end          
    end
end

function [psig, pseps, phprime, pnmu, psn, invK, L, K] = grad_f(sigma, hprime, seps, noise_mu, h, K, Kl, Ksn, heteroseps)
    Ksn_size = size(Ksn);
    nf = Ksn_size(3);
    % calculate K for this round of update
    trueK = K;
    K = K + diag(seps.^2);
    K = preprocess(K, true);
    L = chol(nearestSPD(K))';
    invK = inv_chol(L);
    alpha = invK*(h-noise_mu);
    aat = alpha*alpha';
    ami = aat - invK;
    
    %calculate partial derivative for sigma
    psig = trace(ami*trueK)/sigma;

    %partial derivative for \sigma_\eps
    if heteroseps
        pseps = seps.*diag(ami);
    else 
        pseps = seps*trace(ami);
    end

    %partial derivative for l
    pl = 0.5*trace(ami*Kl);
    phprime = -2*pl/hprime^3;

    %partial derivative for noise_mu
    pnmu = sum((h'-noise_mu)*invK);

    %partial derivative for seps_neuron
    psn = zeros(nf, nf);
%     for i = 1:nf
%         for j = 1:nf
%             psn(i,j) = trace(ami*Ksn(:,:, i,j));
%         end
%     end
    
end
