%% start calculation
truncate = false;
try
    if truncate
        load('result_rel_trunc');
        load('result_mse_trunc');
    else
        load('result_rel');
        load('result_mse');
    end
catch 
    
end
for paramno = 1:4
    if paramno==2
        continue
    end
    nnlist = [200, 500, 1000];
    cd = .1; % connection density
%     nn = 2000; %number of neurons
    mxp = 1; % max # of patterns
    nf = 50; % number of frequency levels
    dim = [nf,nf];
    nent = prod(dim);
    sf = dim(1)/50; %scale factor
    sChange = zeros(dim);

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
    %     sChange = sChange - mean(sChange, 'all');
        folder = "figs_bcm";
    elseif rule == 0 && paramno~=5
        load(['paramset' num2str(paramno)]);
        X = dW1(1:2:end, 1:2:end);
        sChange = X;
        sChange(isnan(sChange)) = 1;
        dim = size(sChange);
        % add structured noise
        for i = 1:nf
            for j = 1:nf
                sChange(i,j) = sqrt(dWvar(i*(100/nf),j*(100/nf)))*randn()+sChange(i,j);
            end
        end
    else
        X = 1/(nn*cd)*.1/meanRate^3*((tmp.^2 - meanRate*tmp)'*...
                                                    (tmp- meanRate));
        sChange = X;
        truncX = X(:, 5:end);
    end

    %% add noise to data
%     noise_index = 1;
%     if noise_index == 0
%         % add uniform noise
%         noise_level = 0.2;
%         sChange = imnoise(sChange, 'gaussian', 0);
%     elseif noise_index == 1
%         % add structured noise
%         for i = 1:nf
%             for j = 1:nf
%                 sChange(i,j) = sqrt(dWvar(i*(100/nf),j*(100/nf)))*randn()+sChange(i,j);
%             end
%         end
%     end 
    
    nature_neuro = false;
    for wc = 1:2
        nn = nnlist(wc);
        np = 1;
        for trial = 1:30
            for nvb = 50:10:50%number of valid observations
                nb = np*nn; %number of observations
                tdW = sChange; %true dW

                %% Network parameters
                k = 2;
                Rate_Novel = gamrnd(k,meanRate/k,nn,np);
                Rate_Novel(Rate_Novel>50) = 50;

                BinConn = (rand(nn)<cd);    % Structural connectivity (it is assumed to be known)
                if nature_neuro
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


                rf = Rate_Novel;
                C = BinConn;


                % sample part of the observations
                % vb stands for valid observation, nvb is the number of valid observations
                vb_index = randperm(nb);
                vb_index = vb_index(1:nvb);
                h = DelW*Rate_Novel + WRec_Novel*Diff_Rate; % + DelW*Diff_Rate;
                h = h(vb_index,:);
                if truncate 
                    ret = ret(:, 5:end);
                    result_rel_trunc(paramno, wc, nvb/10, trial) = norm(ret(:) - truncX(:))/norm(truncX(:));
                    result_mse_trunc(paramno, wc, nvb/10, trial) = immse(ret,truncX);
                    save('result_rel_trunc', 'result_rel_trunc');
                    save('result_mse_trunc', 'result_mse_trunc');
                else 
                    result_rel(paramno, wc, nvb/10, trial) = norm(ret(:) - X(:))/norm(X(:));
                    result_mse(paramno, wc, nvb/10, trial) = immse(ret,X);
                    save('result_rel', 'result_rel');
                    save('result_mse', 'result_mse');
                end
                fprintf("trial = %d nvb = %d\n", trial, nvb);
                norm(ret(:) - X(:))/norm(X(:))
                tmp = WRec_Novel*Diff_Rate;
                mean_diff = abs(mean(tmp) - mean(tmp(vb_index)));
                
                immse(ret,X)
%                   pred = bayesianLinearRegression(nn, np, meanRate, nvb, cd)';
%                   result_reg_rel(wc, nvb/10, trial) = ...
%                       norm(pred(:) - X(:))/norm(X(:));
%                   save('result_reg_rel', 'result_reg_rel');
            end  
            fprintf("wc = %d\n", wc);
        end
    end
end

