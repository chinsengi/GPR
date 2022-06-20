function pred = vanilla_gpr(dim, obi, filledValue)
    nent = prod(dim);
    nn = dim(1);
%     notfilled = setdiff(1:nent, obi);
    filled = length(obi);
    
    % hyperparameter learning
    init_sig = 4*abs(median(filledValue));
    [sigma, l, seps] = vanilla_hyper_infer(dim, obi, filledValue, init_sig); 
    
    K = zeros(filled, filled);
    for i = 1:filled
        for j = 1:filled
            [x1, y1] = unvec(obi(i), nn);
            [x2, y2] = unvec(obi(j), nn);
            K(i,j) = sigma^2*exp(-((x1-x2)^2+(y2-y1)^2)/(2*l^2)); 
        end
    end
    K = K+diag(seps^2*ones(filled,1));
    Ks = zeros(filled, nent);
    for i = 1:filled
        for j = 1:nent
            [x1, y1] = unvec(obi(i), nn);
            [x2, y2] = unvec(j, nn);
            Ks(i,j) = sigma^2*exp(-((x1-x2)^2+(y2-y1)^2)/(2*l^2)); 
        end
    end
%     K = preprocess(K, false);
%     L = chol(nearestSPD(K))';
    alpha = generalizedMinimalResidualMethod(K, filledValue);
    mu_pos =Ks'*alpha;
    pred = zeros(dim);
    for i = 1:nent
        [row, col] = unvec(i, dim(1));
        pred(row, col) = mu_pos(i);
    end
end

