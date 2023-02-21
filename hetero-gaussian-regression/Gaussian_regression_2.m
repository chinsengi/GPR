%% parameter initialization
k = 1; % active sampling rate
dW = zeros(dim); % dW with missing value
filledValue = [];
notfilled = 1:nent;
filled = 0;
obi = []; %observed index

init_sm = 3; 
if init_sm == 1
    % fill in the latin grid
    for i = 1:dim(1)
        dW(i,init(i)) = tdW(i,init(i));
        obi(end+1) = (i-1)*dim(1)+init(i);
        filledValue(end+1) = dW(i,init(i));
        notfilled = setdiff(1:nent, obi);
    end
    filled = filled+dim(1);
elseif init_sm == 2
    % fill in 9 scattered sample
    for i = 1:24:50
        for j = 1:24:50
            dW(i,j) = tdW(i,j);
            obi(end+1) = (j-1)*dim(1)+i;
            filledValue(end+1) = dW(i,j);
            notfilled = setdiff(1:nent, obi);
            filled = filled+1;
        end
    end
else
    % randomly choosing 10 initial sample
    rand_index = randi(50,10,2);
    for index = 1:10
        i = rand_index(index,1);
        j = rand_index(index,2);
        dW(i,j) = tdW(i,j);
        obi(end+1) = (j-1)*dim(1)+i;
        filledValue(end+1) = dW(i,j);
        notfilled = setdiff(1:nent, obi);
        filled = filled+1;
    end
end

%% active sampling and inference
for nitr = 1:105
%     %gradient descent to select model
    if mod(filled,10) == 0 && filled>0
%         dataMean = mean(filledValue);
%         filledValue = filledValue - dataMean;
        dataMean = 0;
        bestlogp = -inf;
        sigma_guess = [0.1, 0.2, 0.05];
        l_guess = [30 40 50];
        seps_guess = [0.05 0.1 0.3];
        for sigma_index = 3:3
            for l_index = 2:2
                for seps_index = 1:1
                    sigma = ones(filled,1)*sigma_guess(sigma_index);
                    l = l_guess(l_index)*ones(filled,1);
                    seps = seps_guess(seps_index)*ones(filled,1); %\sigma_\eps
                    if ~heterosig, sigma = sigma(1); end
                    if ~heterol, l = l(1); end
                    if ~heteroseps, seps = seps(1); end 
                    model_selection_lbgfs;
%                     model_selection_final;
                    if logp>bestlogp
                        bestlogp = logp;
                        sigma_true = sigma;
                        seps_true = seps;
                        l_true = l;
%                         sigma_index
%                         h_index
%                         seps_index
                    end
                end
            end
        end
        l = l_true;
        seps = seps_true;
        sigma = sigma_true;
%         if heteroseps
%             seps = vanilla_gpr(dim, obi, seps); 
%         end
        if heterosig
            sigma = vanilla_gpr(dim, obi, sigma);
        else 
            sigma = ones(nent,1).*sigma(1);
        end
        if heterol
            l = vanilla_gpr(dim, obi, l);
        else
            l = ones(nent,1).*l(1);
        end
    
        % compute the mean value according to conditional gaussian
        K = zeros(filled, filled);
        Ks = zeros(filled, nent);
        for i = 1:filled
            for j = 1:filled
                [x1, y1] = unvec(obi(i), nf);
                [x2, y2] = unvec(obi(j), nf);
                dist = (x1-x2)^2+(y2-y1)^2;
                K(i,j) = sigma(i)*sigma(j)*(2*l(i)*l(j)/(l(i)^2+l(j)^2))...
                                  *exp(-dist/(l(i)^2+l(j)^2));
            end
        end
        for i = 1:filled
            for j = 1:nent
                [x1, y1] = unvec(obi(i), nf);
                [x2, y2] = unvec(j, nf);
                dist = (x1-x2)^2+(y2-y1)^2;
                Ks(i,j) = sigma(i)*sigma(j)*(2*l(i)*l(j)/(l(i)^2+l(j)^2))...
                                  *exp(-dist/(l(i)^2+l(j)^2));        
            end
        end
        if heteroseps
%             seps_filled = seps(obi);
%             seps_filled = seps_filled(:);
            K = K+diag(seps.^2);
        else
            K = K+diag(seps(1)^2*ones(length(K),1));
        end
        alpha = generalizedMinimalResidualMethod(K, filledValue');
%         v = L\Ks;
    %     K_pos = Kss - v'*v; %posterior variance
        mu_pos = Ks'*alpha;
        reconstW = dataMean + reshape(mu_pos, size(dW));
        
        expno = 2;
        hetero = heterosig*4+heterol*2 + heteroseps+1; % write in binary form
        result(hetero, paramno, expno, filled, trial) = norm(reconstW(:)- X(:))/norm(X(:));
        result_logp(hetero,paramno, expno, filled, trial) = bestlogp;
        save('./result/result','result');
        norm(reconstW(:) - X(:))/norm(X(:))
        
        % linear Regression
%         reg_method = 2; % 1 for bayesian regression, 2 for OLS
%         pred4 = SynapticLinearRegression(paramno, obi, filledValue', reg_method, seed, 4);
%         result_rel_reg4_nonoise(paramno, method, filled, trial) = norm(pred4(:) - X(:))/norm(X(:));
%         norm(pred4(:) - X(:))/norm(X(:))
%         pred3 = SynapticLinearRegression(paramno, obi, filledValue', reg_method, seed, 3);
%         result_rel_reg3_nonoise(paramno, method, filled, trial) = norm(pred3(:) - X(:))/norm(X(:));
%         norm(pred3(:) - X(:))/norm(X(:))
%         save('./result/result_rel_reg3_nonoise','result_rel_reg3_nonoise');
%         save('./result/result_rel_reg4_nonoise','result_rel_reg4_nonoise');
        tmp = 1;
    end

    %% active sampling
    if method == 1
        % method 1: varmin
        best_var = inf;
        varmin_result = zeros(length(notfilled));
        for rand_select = 1:length(notfilled)
            pred = rand_select;
            notfilled_pred = notfilled([1:pred-1, pred+1:length(notfilled)]);
            obi_pred = [obi notfilled(pred)];
            Kss = cov(notfilled_pred, notfilled_pred);
            Ks = cov(obi_pred,notfilled_pred);
            K = cov(obi_pred, obi_pred);
            L = chol(K)';
            v = L\Ks;
            K_pos = Kss - v'*v; %posterior variance 
            tr = trace(K_pos);
%             if tr<best_var
%                 selected = notfilled(pred);
%                 best_var = tr;
%             end  
            varmin_result(rand_select) = tr;
        end
        [M,I] = min(varmin_result);
        selected = notfilled(I(1)); % the entry to be queried next
    elseif method == 2
        % method 2: infomax
        best_info = -inf;
        selected = notfilled(1); % the entry to be queried next
        for rand_select = 1:length(notfilled)
            pred = rand_select;
            notfilled_pred = notfilled([1:pred-1, pred+1:length(notfilled)]);
            obi_pred = [obi notfilled(pred)];
            K = cov(obi_pred, obi_pred);
            L = chol(K);
            % not really logdet, discard constant factor
            logdetK = sum(log(diag(L)));
            mi = logdetK;
            if mi>best_info
                selected = notfilled(pred);
                best_info = mi;
            end   
        end
%         infomax_result = zeros(length(notfilled));
%         parfor rand_select = 1:length(notfilled)
%             pred = rand_select;
%             notfilled_pred = notfilled([1:pred-1, pred+1:length(notfilled)]);
%             obi_pred = [obi notfilled(pred)];
%             Kss = cov(notfilled_pred, notfilled_pred);
%             K = cov(obi_pred, obi_pred);
%             L = chol(K);
%             % not really logdet, discard constant factor
%             logdetK = sum(log(diag(L)));
%             L = chol(Kss);
%             logdetKss = sum(log(diag(L)));
%             mi = logdetK+logdetKss;
%             infomax_result(rand_select) = logdetK;   
%         end
%         [M,I] = max(infomax_result);
%         selected = notfilled(I(1));
    elseif method ==3
        %method 3: using ground truth
        best_err = inf;
        selected = 1; % the entry to be queried next
        for rand_select = 1:length(notfilled)
%             pred = ceil(rand()*length(notfilled));
            [x, y] = unvec(notfilled(rand_select), dim(1));
            pred = rand_select;
            filled_pred = filled+1;
            notfilled_pred = notfilled([1:pred-1, pred+1:length(notfilled)]);
            obi_pred = [obi notfilled(pred)];
            filledValue_pred = [filledValue dW(x,y)];
            Kss = cov(notfilled_pred, notfilled_pred);
            Ks = cov(obi_pred,notfilled_pred);
            K = cov(obi_pred, obi_pred);
            L = chol(K)';
            alpha = L'\(L\filledValue_pred(1:filled_pred)');
            mu_pos =Ks'*alpha;
            L = chol(Kss)';
            alpha = L'\(L\mu_pos);
            mu_filled = Ks*alpha;
            reconstW = dW;
            for i = 1:filled_pred
                cur = obi_pred(i);
                [row, col] = unvec(cur, dim(1));
                reconstW(row, col) = mu_filled(i);
            end
            for i = filled_pred+1:nent
                cur = notfilled(i-filled_pred);
                [row, col] = unvec(cur, dim(1));
                reconstW(row, col) = mu_pos(i-filled_pred);
            end
            err = norm(reconstW-X)/2500;
            if err<best_err
                selected = notfilled(pred);
                best_err = err;
            end   
        end
    elseif method == 4
        % picking the entry with largest variance 
        [~, I] = maxk(diag(K_pos), k);
        vg = zeros(dim);
        for i = 1:nent - filled
            [x,y] = unvec(notfilled(i), nf);
            vars = diag(K_pos);
            vg(x,y) = vars(i);
        end
        selected = notfilled(I(k));
    elseif method == 5
        % random sampling
        selected = randi(length(notfilled));
        selected = notfilled(selected);
    end
    if dW(selected)==0
        dW(selected) = tdW(selected);
        obi(end+1) = selected;
        filledValue(end+1) = dW(selected);
        filled = filled + 1;
        notfilled = setdiff(1:nent, obi);
    else
        tmp = 1;
    end  
end