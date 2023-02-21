function pred = SynapticLinearRegression(paramno, obi, h, method, seed, order)
    NumPredictors = 15;
    filled = length(obi);
    [row, col] = ind2sub([50,50], obi); 
    row = row';
    col = col';
    X = zeros(filled, NumPredictors);
    X(:, 1) = 1;
    tmp = [
        row,... %x
        col,... %y
        row.^2,...% x^2
        col.^2,... % y^2
        row.*(col),... % xy
        (row.^2).*(col),... % x^2y
        row.*(col.^2),... %xy^2
        row.*(col.^3),... %xy^3
        (row.^2).*(col.^2),... %x^2y^2
        (row.^3).*(col),... %x^3y
        row.^3,... %x^3
        col.^3,... %y^3
        row.^4,... %x^4
        col.^4,... %y^4
        ];
    X(:,2:end) = tmp;
    
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
            b = regress(h, X);
            MAPBeta = b;
        else
            b = regress(h, X(:, [1:8, 12,13]))';
            MAPBeta = [b(1:8), 0,0,0,b(9:10), 0,0];
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
