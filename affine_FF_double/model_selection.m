% model selection 
function[x_min, logp] = model_selection(X0, model_selection_method, h,...
                                vb_index, C, rf, rfpreI, rfpreE, postgrad)
    optmethod = 1;
    checkGradient = false;
    if model_selection_method == 1
        grad_f = @grad_f_ML;
    elseif model_selection_method == 2
        optmethod = 3;
    else
        grad_f = @grad_f_cval;
    end


    loss = @(X)grad_f(X, h, vb_index, C, rf, rfpreI, rfpreE);
    if nargin>9
        loss = @(X)grad_f(X, h, vb_index, C, rf, rfpreI, rfpreE, postgrad);
    end

    if optmethod ==1
        options = optimoptions('fmincon','Algorithm','sqp',...
            'SpecifyObjectiveGradient',true,...
            'Display', 'off',...
            'CheckGradients', checkGradient,...
            'FiniteDifferenceType', 'central',...
            'MaxIterations', 25 );
        lb = -10*ones(size(X0)); %lower bound 
        lb(2) = 2; % the length_scale cannot be too small
        lb(5) = 2; % the length_scale cannot be too small
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
        logp = -grad_f_ML(x_min, h, vb_index, C, rf, rfpreI, rfpreE, postgrad);
    else
        logp = -grad_f_ML(x_min, h, vb_index, C, rf, rfpreI, rfpreE);
    end
end

%gradient function for cross validation
function [f, g] = grad_f_cval(X, h, vb_index, C, rf, rfpreI, rfpreE, postgrad)
    input_noise = nargin > 7;
    % calculate K for this round of update
    if input_noise
        [K, KlE, KlI, KsnE, KsnI, KsigE, KsigI, gtg] =...
                   Kgen(X, vb_index, C, rf, rfpreI, rfpreE,postgrad);
    else
        [K, KlE, KlI, KsnE, KsnI, KsigE, KsigI] =...
                   Kgen(X, vb_index, C, rf, rfpreI, rfpreE);
    end
    
    try
        L = chol(K)';
    catch
        L = chol(nearestSPD(K))';
    end
    invK = inv_chol(L);
    alpha = invK*(h-mu);
    dinvK = diag(invK);
    diff = alpha./dinvK-mu;
    f = -0.5*sum(log(dinvK))+ sum(diff.^2.*dinvK)/2;  
%     f = -0.5*sum(log(dinvK))+ sum(diff.^2.*dinvK)/2;
    if ( nargout > 1 )
        %intermediate function for partial derivative calculation
        derivTmp = @(ztheta)...
            (alpha-mu*dinvK).*(ztheta*alpha)...
            - 0.5*(1+(alpha.^2-(mu*dinvK).^2)./dinvK).*diag(ztheta*invK);
        
        %calculate partial derivative for sigma
        zsig = invK*Ksig;
        psig = derivTmp(zsig);
        psig = sum(psig./dinvK);
        
        %partial derivative for \sigma_\eps
        zseps = 2*seps*invK;
        pseps = derivTmp(zseps);
        pseps = sum(pseps./dinvK);


        %partial derivative for l
        zl = invK*Kl;
        pl = derivTmp(zl);
        pl = sum(pl./dinvK);
        
        %partial derivative for mu
        pnmu = sum((alpha - mu*dinvK).*(1+(sum(invK, 2))./dinvK));

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
function [f, g] = grad_f_ML(X, h, vb_index, C, rf, rfpreI, rfpreE, postgrad)
    input_noise = nargin > 7;
    nvb = numel(vb_index);
    [sigmaE, lE, snE, sigmaI, lI, snI, innoise] =...
                                                     extract_param(X);
    if input_noise
        [K, KlE, KlI, KsnE, KsnI, KsigE, KsigI, gtg] =...
                         Kgen(X, vb_index, C, rf, rfpreI, rfpreE,postgrad);
    else
        [K, KlE, KlI, KsnE, KsnI, KsigE, KsigI] =...
                                 Kgen(X, vb_index, C, rf, rfpreI, rfpreE);
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
        invK = inv_chol(L);
        alpha = K\(h);
        aat = alpha*alpha';
        ami = aat - invK/theta;

        %calculate partial derivative for sigma
        psigE = 0.5*trace(ami*KsigE);
        psigI = 0.5*trace(ami*KsigI);

        %partial derivative for l
        plE = 0.5*trace(ami*KlE);
        plI = 0.5*trace(ami*KlI);
        
        %partial derivative for seps_neuron
        psnE = 0.5*trace(ami*KsnE);
        psnI = 0.5*trace(ami*KsnI);
        
        %partial derivative for input noise variance
        if input_noise
            pinnoise = innoise*trace(ami*gtg);
            g = [psigE; plE; psnE; psigI; plI; psnI; pinnoise];
        else
%             g = [psig; pl; pseps; pnmu; psn];
            g = [psigE; plE; psnE; psigI; plI; psnI];
        end
        g = -g.*exp(X);
    end
end