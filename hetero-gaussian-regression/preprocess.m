function K = preprocess(K, fast) %deal with numerical issue
    K(isnan(K)) = 0;
    K(K == inf) = 0;
    correction = 1e-7;
    if ~fast
        while cond(K)>1e6
            K = K+diag(correction*ones(length(K),1));
            correction = correction*1.5;
        end
    else
        dK = decomposition(K);
        if isIllConditioned(dK)
            K = K+diag(ones(length(K),1));
        end
    end
end