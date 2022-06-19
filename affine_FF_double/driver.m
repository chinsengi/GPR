%% load result
% load('result_nmseE');
% load('result_relerrE');
% load('result_nmseI');
% load('result_relerrI');
% load('save_seed');

%% start calculation
cd = .1; % connection density
mxp = 1; % max # of patterns
nf = 50; % number of frequency levels 
dim = [nf,nf];
nent = prod(dim);       
sChange = zeros(dim);
input_noise = false; % whether input noise is taken into consideration

for paramno = 1:1
    if paramno==2
        continue
    end

    rule = 0; % 0 stands for calcium-based rule, 1 stands for BCM based rule
    if rule == 1
        theta1 = 0;
        theta2 = 8; %BCM rule
        for i = 1:nn
            g = i/sf - 25;
            for j = 1:nn
                if i <= theta1
                    f = 0;
                else
                    f = (j/sf-theta1)*(j/sf-theta2);
                end
                sChange(i,j) = f*g;
            end
        end
        sChange = sChange/10000;
        folder = "figs_bcm";
    elseif rule == 0 && paramno ~= 0
%         load(['paramset' num2str(paramno)]);
        load('param1max70.mat');
%         X = dW1(1:2:end, 1:2:end);
        X = dW1;
        sChange = X;
        sChange(isnan(sChange)) = 1;
        dim = size(sChange);
    end

    if paramno ==0 
        heteroseps = false;
        tmp = 1:50;
        X = 1/(nn*cd)*.1/meanRate^3*((tmp.^2 - meanRate*tmp)'*...
                                                        (tmp- meanRate));
        sChange = X;
        truncX = X(:, 5:end);
    else
        heteroseps = false;
    end

    nn = 500;
    np = 1;
    Xtrue = X;
    X = 1/(nn*cd)*(X-0.5)/.5;
    expno = 1;
    % exp 1: ordinary double FF (same map norm, 4:1 EI ratio)
    % exp 2: EI network(2:1 EI map norm ratio, 1:1 EI ratio, post dist uniform)
    % exp 3: EI network(1:2 EI map norm ratio)
    % exp 4: EI network(2:3 EI map norm ratio)
    % exp 5: EI network(3:2 EI map norm ratio)
    % exp 6: EE network(same map norm)
    % exp 7: EI network(4:1 EI map norm ratio)
    % exp 8: EI network(3:1 EI map norm ratio)
    % exp 9: EI network(1:1 EI map norm ratio)
    % exp 10: EI network(1:3 EI map norm ratio)
    % exp 11: EI network(1:4 EI map norm ratio)
    for expno = 1:1
        for trial = 1:50
            %set seed
            rng shuffle
            seed = randi(100000);
            save_seed(expno, paramno+1, trial) = seed;
            save('save_seed', 'save_seed');   
%             seed = 38145 ;
            rng(seed);
            % add structured noise
            for i = 1:length(dWvar)
                for j = 1:length(dWvar)
                    sChange(i,j) = sqrt(dWvar(i,j))*randn() + Xtrue(i,j);
    %                     sChange(i,j) = Xtrue(i,j)+randn()*0.04;
                end
            end
            sChange = 1/(nn*cd)*(sChange - 0.5)/.5;
            for nvb = 10:10:100 %number of valid observations
                for model_selection_method = 1:1
                    infer_rank1; 
                    % ir = inference range
                    % set inference range
%                     medianE = round(median(pdE));
%                     medianI = round(median(pdI));
%                     irEx = max(1, (medianE(1) - 15)):min(50, (medianE(1)+15));
%                     irEy = max(1, (medianE(2) - 15)):min(50, (medianE(2)+15));
%                     irIx = max(1, (medianI(1) - 15)):min(50, (medianI(1)+15));
%                     irIy = max(1, (medianI(2) - 15)):min(50, (medianI(2)+15));
                    irEx = 5:40; irEy = 5:40;
                    if expno>1
                        irIx = 5:40; irIy = 5:40;
                    else
                        irIx = 10:45; irIy = 10:45;
                    end
                    trueE = X(irEx, irEy);
                    trueI = XI(irIx, irIy);
                    predE = dWE(irEx, irEy);
                    predI = dWI(irIx, irIy);
                    result_relerrE(expno, paramno+1, model_selection_method,...
                                                             nvb/10, trial)...
                           = relativeError(predE, trueE);
                    result_nmseE(expno, paramno+1, model_selection_method,...
                                                             nvb/10, trial)...
                           = nmse(predE, trueE);
                    result_relerrI(expno, paramno+1, model_selection_method,...
                                                             nvb/10, trial)...
                           = relativeError(predI, trueI);
                    result_nmseI(expno, paramno+1, model_selection_method,...
                                                             nvb/10, trial)...
                           = nmse(predI, trueI);
    %                     if model_selection_method == 1
    %                         paramtrained(trial, nvb/10, :) = param_true;
    %                     end
                    recovered_I(expno, paramno+1, model_selection_method,...
                                                     nvb/10, trial,:,:) = dWI;
                    recovered_E(expno, paramno+1, model_selection_method,...
                                                     nvb/10, trial,:,:) = dWE;
                    relativeError(predE, trueE)
                    relativeError(predI, trueI)
                    tmp = 1;
    %                 nmse(predE, trueE)
    %                 nmse(predI, trueI)
                end
%                 save('result_relerrE', 'result_relerrE');
%                 save('result_nmseE', 'result_nmseE');
%                 save('result_relerrI', 'result_relerrI');
%                 save('result_nmseI', 'result_nmseI');
    
                %linear regression performance for comparison
    %             reg_method = 2; % 1 for bayesian, 2 for ordinary MLS
    %             pred4 = bayesianLinearRegression(nn, np, meanRate, nvb, cd, sChange, paramno, reg_method, seed, 4);
    %             result_reg_rel4(paramno+1, wc, nvb/10, trial) = ...
    %                 norm(pred4(:) - X(:))/norm(X(:));
    %             norm(pred4(:) - X(:))/norm(X(:))
    %             pred3 = bayesianLinearRegression(nn, np, meanRate, nvb, cd, sChange, paramno, reg_method, seed, 3);
    %             result_reg_rel3(paramno+1, wc, nvb/10, trial) = ...
    %                 norm(pred3(:) - X(:))/norm(X(:));
    %             norm(pred3(:) - X(:))/norm(X(:))
    %             save('result_reg_rel3', 'result_reg_rel3');
    %             save('result_reg_rel4', 'result_reg_rel4');
                if nvb == 100
                    tmp = 1;
                end
                fprintf("expno = %d, paramno = %d, trial = %d nvb = %d\n", expno, paramno, trial, nvb);
            end  
%             save('result_relerrE', 'result_relerrE');
%             save('result_nmseE', 'result_nmseE');
%             save('result_relerrI', 'result_relerrI');
%             save('result_nmseI', 'result_nmseI');
        end
    end
end

function err = relativeError(A, B)
    err = norm(A(:) - B(:))/norm(B(:));
end

function rsq = rSquared(A, B)
    rsq = 1 - norm(A(:) - B(:))^2/norm(B(:) - mean(B(:)))^2;
end

function NMSE = nmse(A, B)
    NMSE = norm(A(:) - B(:))^2/norm(B(:) - mean(B(:)))^2;
end