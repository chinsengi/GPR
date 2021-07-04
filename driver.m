%% start calculation
truncate = false;
use_gpr = true;
use_reg = true;
try
    if use_gpr
        load('result_rel');
        load('result_mse');
    end
    if use_reg
        load('result_reg_rel');
    end
    load('mean_diff');
catch 
    
end
nnlist = [200, 500, 1000];
cd = .1; % connection density
mxp = 1; % max # of patterns
nf = 50; % number of frequency levels
meanRate = 20;  % Average firing rate before learning
dim = [nf,nf];
nent = prod(dim);
sf = dim(1)/50; %scale factor
sChange = zeros(dim);
input_noise = false; % whether input noise is taken into consideration

for paramno = 1:1
    for wc = 2:2
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
            load(['paramset' num2str(paramno)]);
            X = dW1(1:2:end, 1:2:end);
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
        
        nn = nnlist(wc);
        np = 1;
        Xtrue = X;
        X = 1/(nn*cd)*(X-0.5)/.5;
        for trial = 1:20
            %set seed
            seed = randi(100000);
%             save_seed(paramno+1, trial, wc) = seed;
%             save('save_seed', 'save_seed');
%             seed = save_seed(paramno+1, trial, wc);
%             seed = 72155;
%             rng(seed);
            % add structured noise
            for i = 1:nf
                for j = 1:nf
%                     sChange(i,j) = sqrt(dWvar(i*(100/nf),j*(100/nf)))*randn()...
%                                                                +Xtrue(i,j);
                    sChange(i,j) = Xtrue(i,j);
                end
            end
            sChange = 1/(nn*cd)*(sChange - 0.5)/.5;
            for nvb = 10:10:100 %number of valid observations
%                 for model_selection_method = 1:3
                    infer_rank1;  
                    if truncate 
                        ret = ret(:, 5:end);
                        result_rel_trunc(paramno+1, wc, nvb/10, trial) = norm(ret(:) - truncX(:))/norm(truncX(:));
                        result_mse_trunc(paramno+1, wc, nvb/10, trial) = immse(ret,truncX);
                        save('result_rel_trunc', 'result_rel_trunc');
                        save('result_mse_trunc', 'result_mse_trunc');
                    else 
                        result_rel_fam(paramno+1, model_selection_method, nvb/10, trial) = norm(ret(:) - X(:))/norm(X(:));
                        result_logp(paramno+1, model_selection_method, nvb/10, trial) = bestlogp;
                        result_mse(paramno+1, wc, nvb/10, trial) = immse(ret,X);
%                         save('result_rel_fam', 'result_rel_fam');
%                         save('result_logp', 'result_logp');
    %                     save('result_rel_input_vanilla', 'result_rel_input_vanilla');
    %                     save('result_rel', 'result_rel');
    %                     save('result_mse', 'result_mse');
                    end
                    norm(ret(:) - X(:))/norm(X(:))
%                 end
                %linear regression performance for comparison
%                     reg_method = 2; % 1 for bayesian, 2 for ordinary MLS
%                     pred4 = bayesianLinearRegression(nn, np, meanRate, nvb, cd, sChange, paramno, reg_method, seed, 4);
%                     result_reg_rel4(paramno+1, wc, nvb/10, trial) = ...
%                         norm(pred4(:) - X(:))/norm(X(:));
%                     norm(pred4(:) - X(:))/norm(X(:))
%                     pred3 = bayesianLinearRegression(nn, np, meanRate, nvb, cd, sChange, paramno, reg_method, seed, 3);
%                     result_reg_rel3(paramno+1, wc, nvb/10, trial) = ...
%                         norm(pred3(:) - X(:))/norm(X(:));
%                     norm(pred3(:) - X(:))/norm(X(:))
%                     save('result_reg_rel3', 'result_reg_rel3');
%                     save('result_reg_rel4', 'result_reg_rel4');
                fprintf("paramno = %d, trial = %d nvb = %d\n", paramno, trial, nvb);
            end  
        end
    end
end
