%% load result
% load('result_nmse');
% load('result_relerr');
% load('result_logp');
% load('result_reg_rel3');
% load('result_reg_rel4');
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

    if paramno ~= 0
%         load(['paramset' num2str(paramno)]);
%         X = dW1(1:2:end, 1:2:end);
        load('param1max70.mat');
%         X = dW1(1:2:end, 1:2:end);
        X = dW1;
        sChange = X;
        sChange(isnan(sChange)) = 1;
        dim = size(sChange);
    elseif paramno ==0 
        heteroseps = false;
        tmp = 1:50;
        X = 1/(nn*cd)*.1/meanRate^3*((tmp.^2 - meanRate*tmp)'*...
                                                        (tmp- meanRate));
        sChange = X;
        truncX = X(:, 5:end);
    end

    np = 1;
    Xtrue = X;
    X = 1/50*(X-0.5)/.5;
    
    % expno = 
    % 1: normal inference
    % 2: missing connection 10%
    % 3: input noise w/o treatment, noise_level = 0.1
    % 4: input noise with treatment, noise_level = 0.1
    % 5: input noise w/o treatment, noise_level = 0.5
    % 6: input noise with treatment, noise_level = 0.5
    % 7: missing connection 20%
    % 8: active protocal
    % 9: 100 pre-synaptic neurons
    % 10: 200 pre-synaptic neurons
    % 11: 300 pre-synaptic neurons
    % 12: 400 pre-synaptic neurons
    for expno = 1:1
        if expno== 8 
            continue 
        end

        if expno<9
            nn = 500;
        elseif expno<=12
            nn = (expno - 8)*100;
        end
        for trial = 1:50
            %set seed
            rng shuffle
            seed = randi(100000);
            save_seed(expno, paramno+1, trial) = seed;
%             save('save_seed', 'save_seed');
%             seed = 32082;
            rng(seed);
            % add structured noise
            for i = 1:70
                for j = 1:70
%                     sChange(i,j) = sqrt(dWvar(i*(100/nf),j*(100/nf)))*randn()...
%                                                                +Xtrue(i,j);
                      sChange(i,j) = sqrt(dWvar(i,j))*randn() + Xtrue(i,j);
%                     sChange(i,j) = Xtrue(i,j)+randn()*0.04;
                end
            end
            sChange = 1/50*(sChange - 0.5)/.5;
            
            input_noise = (expno==4 || expno==6);
            
            for nvb = 10:10:100 %number of valid observations
                irx = 5:40; iry = 5:40; % ir = inference range
                % model_selection_method is 1 -> log likelihood
                % 2 -> no hyperparameter training
                % 3 -> cross validation
                for model_selection_method = 1:1
                    infer_rank1;  
                    Xp = X(irx, iry); ret = ret(irx, iry);
                    result_relerr(expno, paramno+1, model_selection_method,...
                        nvb/5, trial) = relativeError(ret, Xp);
                    result_logp(expno, paramno+1, model_selection_method,...
                        nvb/5, trial) = bestlogp;
                    relativeError(ret, Xp)
                end
                %linear regression performance for comparison
                reg_method = 2; % 1 for bayesian, 2 for ordinary MLS
%                 pred4 = bayesianLinearRegression(nn, np, npost, meanRate, nvb, cd, sChange, reg_method, seed, 4, true);
%                 pred4 = pred4(irx, iry);
%                 result_reg_rel4(paramno+1, wc, nvb/10, trial) = ...
%                     norm(pred4(:) - X(:))/norm(X(:));
                pred3 = bayesianLinearRegression(nn, np, npost, meanRate, nvb, cd, sChange, reg_method, seed, 3, true);
                pred3 = pred3(irx, iry);
%                 result_reg_rel3(paramno+1, wc, nvb/10, trial) = ...
%                     relativeError(pred3, Xp);
                relativeError(pred3, Xp)
%                 save('result_reg_rel3', 'result_reg_rel3');
%                 save('result_reg_rel4', 'result_reg_rel4');
%                 save('result_relerr', 'result_relerr');
%                 save('result_nmse', 'result_nmse');
                fprintf("exp: %d, paramno = %d, trial = %d nvb = %d\n", expno, paramno, trial, nvb);
            end  
        end
    end
end

function err = relativeError(A, B)
    err = norm(A(:) - B(:))/ norm(B(:));
end

function rsq = rSquared(A, B)
    rsq = 1-norm(A(:) - B(:))^2/norm(B(:) - mean(B(:)))^2;
end

function NMSE = nmse(A, B)
    NMSE = norm(A(:) - B(:))^2/norm(B(:) - mean(B(:)))^2;
end