% model selection 
function[x_min, logp] = model_selection(X0, model_selection_method, h,...
                                vb_index, C, rf, rfpre, postgrad)
    optmethod = 1;
    checkGradient = false;
    if model_selection_method == 1
        grad_f = @grad_f_ML;
    elseif model_selection_method == 2
        optmethod = 3;
    else
        grad_f = @grad_f_cval;
    end


    loss = @(X)grad_f(X, h, vb_index, C, rf, rfpre);
    if nargin>7
        loss = @(X)grad_f(X, h, vb_index, C, rf, rfpre, postgrad);
    end

    if optmethod ==1
        options = optimoptions('fmincon','Algorithm','interior-point',...
            'SpecifyObjectiveGradient',true,...
            'Display', 'off',...
            'CheckGradients', checkGradient,...
            'FiniteDifferenceType', 'central',...
            'MaxIterations', 40 );
        lb = -10*ones(size(X0)); %lower bound 
        lb(2) = 2; % the length_scale cannot be too small
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
    if nargin>8
        logp = -grad_f_ML(x_min, h, vb_index, C, rf, rfpre, postgrad);
    else
        logp = -grad_f_ML(x_min, h, vb_index, C, rf, rfpre);
    end
end

%gradient function for cross validation
function [f, g] = grad_f_cval(X, h, vb_index, C, rf, rfpre, postgrad)
    input_noise = nargin > 6;
    [~, ~,  ~, innoise] = extract_param(X);
    if input_noise
        [K, Kl, Ksn, Ksig, gtg] = Kgen(X, vb_index, C, rf, rfpre,...
                                                         postgrad);
    else
        [K, Kl, Ksn, Ksig] = Kgen(X, vb_index, C, rf, rfpre);
    end
    try
        L = chol(K)';
    catch
        L = chol(nearestSPD(K))';
    end
    invK = inv_chol(L);
    alpha = invK*h;
    dinvK = diag(invK);
    diff = alpha./dinvK;
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
%             g = [psig; pl; pseps; pnmu; psn];
            g = [psig; pl; pseps; pnmu; psn];
        end
        g = -g.*exp(X);
    end
end

%gradient function for MAP
function [f, g] = grad_f_ML(X, h, vb_index, C, rf, rfpre, postgrad)
    input_noise = nargin > 6;
    [sig, l, sn, innoise] = extract_param(X);
    if input_noise
        [K, Kl, Ksn, Ksig, gtg] = Kgen(X, vb_index, C, rf, rfpre,...
                                                         postgrad);
    else
        [K, Kl, Ksn, Ksig] = Kgen(X, vb_index, C, rf, rfpre);
    end
    % calculate K for this round of update
    if any(isnan(K(:))) || any(isinf(K(:)))
        K(isnan(K)) = 0;
        K(isinf(K)) = 0;
    end
    try
        L = chol(K)';
    catch
        L = chol(nearestSPD(K))';
    end
    beta = L\h;
    theta = 1;
    f = 0.5*(beta'*beta)+sum(log(diag(L)))/theta;
    
    if ( nargout > 1 )
        invK = inv(L)'*inv(L);
%         alpha = K\h;
        alpha = invK*h;
        aat = alpha*alpha';
        ami = aat - invK/theta;

        %calculate partial derivative for sigma
        psig = 0.5*trace(ami*Ksig);

        %partial derivative for l
        pl = 0.5*trace(ami*Kl);

        %partial derivative for seps_neuron
        psn = 0.5*trace(ami*Ksn);
        
        %partial derivative for input noise variance
        if input_noise
            pinnoise = innoise*trace(ami*gtg);
%             pinnoise = 0;
            g = [psig; pl; psn; pinnoise];
        else
            g = [psig; pl; psn];
%             g = [psig; pl; 0; 0; psn];
        end
        g = -g.*exp(X);
    end
end