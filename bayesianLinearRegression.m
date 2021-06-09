function pred = bayesianLinearRegression(nn, np, meanRate, nvb, cd, sChange, paramno, method, seed, order)
    rng(seed);

    nb = np*nn; %number of observations

    %% Network parameters
    k = 2;
    Rate_Novel = gamrnd(k,meanRate/k,nn,np);
    Rate_Novel(Rate_Novel>50) = 50;

    BinConn = (rand(nn)<cd);    % Structural connectivity (it is assumed to be known)
    if paramno==0
        StrengthConn = 1/(nn*cd)*(2*rand(nn)+4);   % Strength of connection before learning
        WRec_Novel = 1/meanRate*BinConn.*StrengthConn;

        IExt = Rate_Novel-WRec_Novel*Rate_Novel;            % External input (assumed to be unchanged with learning)
        %% Synaptic plasticity and firing rate after learning for each pattern
        DelW_Strength = zeros(nn,nn,np);
        DelW_Strength_total = zeros(nn);
        for i = 1:np
            DelW_Strength(:,:,i) = BinConn.*(1/(nn*cd)*.1/meanRate^3*...
                (Rate_Novel(:,i).^2-meanRate*Rate_Novel(:,i))*(Rate_Novel(:,i)'-meanRate));
            DelW_Strength_total = DelW_Strength_total+DelW_Strength(:,:,i);
        end

        DelW = DelW_Strength_total;
        % Diff_Rate = (eye(nn) - WRec_Novel - DelW)\(DelW*Rate_Novel);
        WRec_Fam = WRec_Novel + DelW;
        Rate_Fam = (eye(nn) - WRec_Fam)\IExt;
        %%  Inference on the post-synaptic dependence
        Diff_Rate = Rate_Fam-Rate_Novel;
    else
        Rate_Novel = ceil(Rate_Novel);
        StrengthConn = 1/(nn*cd)*(2*rand(nn)+4);   % Strength of connection before learning
        WRec_Novel = 1/meanRate*BinConn.*StrengthConn;

        IExt = Rate_Novel-WRec_Novel*Rate_Novel;            % External input (assumed to be unchanged with learning)
        DelW = zeros(nn);
        for i = 1:nn
            for j = 1:nn
                if BinConn(i,j)~=0
                    DelW(i,j) = sChange(Rate_Novel(i), Rate_Novel(j));
                end
            end
        end 
        WRec_Fam = WRec_Novel + DelW;
        Rate_Fam = (eye(nn) - WRec_Fam)\IExt;
        %%  Inference on the post-synaptic dependence
        Diff_Rate = Rate_Fam-Rate_Novel;
    end

    NumPredictors = 16;
    % sample part of the observations
    % vb stands for valid observation, nvb is the number of it
    vb_index = randperm(nb);
    vb_index = vb_index(1:nvb);
    h = DelW*Rate_Novel + WRec_Novel*Diff_Rate; % + DelW*Diff_Rate;
    h = h(vb_index,:);
    X = zeros(nvb, NumPredictors);
    X(:, 1) = 1;
    for i = 1:nvb
        ind = vb_index(i);
        tmp = [Rate_Novel,... %constant
            Rate_Novel(ind)*Rate_Novel,... %x
            Rate_Novel.^2,... %y
            Rate_Novel(ind)^2*Rate_Novel,...% x^2
            Rate_Novel.^3,... % y^2
            Rate_Novel(ind)*(Rate_Novel.^2),... % xy
            Rate_Novel(ind)^2*(Rate_Novel.^2),... % x^2y
            Rate_Novel(ind)*(Rate_Novel.^3),... %xy^2
            Rate_Novel(ind)*(Rate_Novel.^4),... %xy^3
            Rate_Novel(ind)^2*(Rate_Novel.^3),... %x^2y^2
            Rate_Novel(ind)^3*(Rate_Novel.^2),... %x^3y
            Rate_Novel(ind)^3*Rate_Novel,... %x^3
            Rate_Novel.^4,... %y^3
            Rate_Novel(ind)^4*Rate_Novel,... %x^4
            Rate_Novel.^5,... %y^4
            ];
        X(i,2:end) = BinConn(ind,:)*tmp;
    end
    
    noise_mu = 1;
    if method == 1
        X = X(:, 3:end);
        InterceptPriorMean = 0;
        InterceptPriorSigma = 1e-3;
%         BetaPriorMean = [0, 0, 1/(nn*cd)*.1/meanRate,...
%         0,-1/(nn*cd)*.1/meanRate^2, -1/(nn*cd)*.1/meanRate^2,...
%         0,1/(nn*cd)*.1/meanRate^3]';
        BetaPriorMean  = 1e-6;
        BetaPriorSigma = 1e-6;
        LogNoiseVarianceMean = 0;
        LogNoiseVarianceSigma = 1e-3;
        logpdf = @(Parameters)logPosterior(Parameters,X,h, ...
            InterceptPriorMean,InterceptPriorSigma, ...
            BetaPriorMean,BetaPriorSigma, ...
            LogNoiseVarianceMean,LogNoiseVarianceSigma);
        Intercept = 0;
        Beta = randn(7,1)*1e-6;
        LogNoiseVariance = 0.02;
        startpoint = [Intercept;Beta;LogNoiseVariance];
        smp = hmcSampler(logpdf,startpoint,'NumSteps',50);
        [MAPpars,fitInfo] = estimateMAP(smp,'VerbosityLevel',0);
        MAPIntercept = MAPpars(1);
        MAPBeta = MAPpars(1:end-1);
        MAPLogNoiseVariance = MAPpars(end);
    else 
%         b = ridge(h,X(:, 2:end), [0.001, 0.01, 0.1, 0.2, 0.5]);
%         MAPBeta = b(1:end, 5);
        if order == 4
            b = regress(h-noise_mu, X(:,2:end));
            MAPBeta = b;
        else
            b = regress(h-noise_mu, X(:,[2:9,13, 14]))';
            MAPBeta = [b(1:8), 0,0,0,b(9), b(10),0,0];
        end
    end
    cpre = (1:50)';
    cpost = (1:50);
    pred = MAPBeta(1) + ones(50,1)*cpost*MAPBeta(2)+cpre*ones(1,50)*MAPBeta(3)...
            +ones(50,1)*(cpost.^2)*MAPBeta(4) + (cpre.^2)*ones(1,50)*MAPBeta(5)...
            +cpre*cpost*MAPBeta(6)+cpre*(cpost.^2)*MAPBeta(7)+(cpre.^2)*cpost*MAPBeta(8)...
            +(cpre.^3)*cpost*MAPBeta(9)+(cpre.^2)*(cpost.^2)*MAPBeta(10)+cpre*(cpost.^3)*MAPBeta(11)...
            +ones(50,1)*cpost.^3*MAPBeta(12)+cpre.^3*ones(1,50)*MAPBeta(13)+...
            ones(50,1)*cpost.^4*MAPBeta(14)+cpre.^4*ones(1,50)*MAPBeta(15);
    pred = pred';
%     figure(1);
%     surf(pred);
%     
%     MAPBeta = [0, 0,1/(nn*cd)*.1/meanRate,...
%         0,-1/(nn*cd)*.1/meanRate^2, -1/(nn*cd)*.1/meanRate^2,...
%         0,1/(nn*cd)*.1/meanRate^3];
%     truth =  MAPBeta(1) + ones(50,1)*cpost*MAPBeta(2)+cpre*ones(1,50)*MAPBeta(3)...
%             +ones(50,1)*(cpost.^2)*MAPBeta(4) + (cpre.^2)*ones(1,50)*MAPBeta(5)...
%             +cpre*cpost*MAPBeta(6)+cpre*(cpost.^2)*MAPBeta(7)+(cpre.^2)*cpost*MAPBeta(8);
% 
%     figure(2);
%     surf(truth);
%     MAPBeta = [0, 1/(nn*cd)*.1/meanRate,0,...
%         -1/(nn*cd)*.1/meanRate^2,0,-1/(nn*cd)*.1/meanRate^2,...
%         1/(nn*cd)*.1/meanRate^3,0]; 
%     MAPBeta = [0, 0, 1/(nn*cd)*.1/meanRate,...
%         0,-1/(nn*cd)*.1/meanRate^2, -1/(nn*cd)*.1/meanRate^2,...
%         0,1/(nn*cd)*.1/meanRate^3];
%     norm(h - X*MAPBeta(2:end)')
%     norm(h - X*b)
    
%     trueh = h;
%     MAPBeta = [0, 1/(nn*cd)*.1/meanRate,0,...
%         -1/(nn*cd)*.1/meanRate^2,0, -1/(nn*cd)*.1/meanRate^2,...
%         1/(nn*cd)*.1/meanRate^3,0];
%     cpre = Rate_Novel;
%     cpost = Rate_Novel';
%     DelW_Strength = zeros(nn,nn,np);
%     DelW_Strength_total = zeros(nn);
%     for i = 1:np
%         DelW_Strength(:,:,i) = BinConn.*(ones(nn,1)*cpost*MAPBeta(3)+cpre*ones(1,nn)*MAPBeta(2)...
%                 +ones(nn,1)*(cpost.^2)*MAPBeta(5) + (cpre.^2)*ones(1,nn)*MAPBeta(4)...
%                 +cpre*cpost*MAPBeta(6)+cpre*(cpost.^2)*MAPBeta(8)+(cpre.^2)*cpost*MAPBeta(7));
% 
%         DelW_Strength_total = DelW_Strength_total+DelW_Strength(:,:,i);
%     end
% 
%     DelW = DelW_Strength_total;
%     h = DelW*Rate_Novel;% + WRec_Novel*Diff_Rate; % + DelW*Diff_Rate;
%     h = h(vb_index,:);
%     mean(trueh - h)
end  

%% Functions for Computing Posterior Distribution
% The |logPosterior| function returns the logarithm of the product of a
% normal likelihood and a normal prior for the linear model. The input
% argument |Parameter| has the format |[Intercept;Beta;LogNoiseVariance]|.
% |X| and |Y| contain the values of the predictors and response,
% respectively.
%
% The |normalPrior| function returns the logarithm of the multivariate
% normal probability density with means |Mu| and standard deviations
% |Sigma|, specified as scalars or columns vectors the same length as |P|.
% The second output argument is the corresponding gradient.
function [logpdf, gradlogpdf] = logPosterior(Parameters,X,Y, ...
    InterceptPriorMean,InterceptPriorSigma, ...
    BetaPriorMean,BetaPriorSigma, ...
    LogNoiseVarianceMean,LogNoiseVarianceSigma)


    % Unpack the parameter vector
    Intercept        = Parameters(1);
    Beta             = Parameters(2:end-1);
    LogNoiseVariance = Parameters(end);
    % Compute the log likelihood and its gradient
    Sigma                   = sqrt(exp(LogNoiseVariance));
    Mu                      = X*Beta + Intercept;
    Z                       = (Y - Mu)/Sigma;
    loglik                  = sum(-log(Sigma) - .5*log(2*pi) - .5*Z.^2);
    gradIntercept1          = sum(Z/Sigma);
    gradBeta1               = X'*Z/Sigma;
    gradLogNoiseVariance1	= sum(-.5 + .5*(Z.^2));
    % Compute log priors and gradients
    [LPIntercept, gradIntercept2]           = normalPrior(Intercept,InterceptPriorMean,InterceptPriorSigma);
    [LPBeta, gradBeta2]                     = normalPrior(Beta,BetaPriorMean,BetaPriorSigma);
    [LPLogNoiseVar, gradLogNoiseVariance2]  = normalPrior(LogNoiseVariance,LogNoiseVarianceMean,LogNoiseVarianceSigma);
    logprior                                = LPIntercept + LPBeta + LPLogNoiseVar;
    % Return the log posterior and its gradient
    logpdf               = loglik + logprior;
    gradIntercept        = gradIntercept1 + gradIntercept2;
    gradBeta             = gradBeta1 + gradBeta2;
    gradLogNoiseVariance = gradLogNoiseVariance1 + gradLogNoiseVariance2;
    gradlogpdf           = [gradIntercept;gradBeta;gradLogNoiseVariance];
end

function [logpdf,gradlogpdf] = normalPrior(P,Mu,Sigma)
Z          = (P - Mu)./Sigma;
logpdf     = sum(-log(Sigma) - .5*log(2*pi) - .5*(Z.^2));
gradlogpdf = -Z./Sigma;
end
