%% initialize parameter
sigma = 1;
l = 0.0007;
h = sqrt(1/l);
seps = reshape(sqrt(dWvar(1:2:100, 1:2:100)'),nent,1); %\sigma_\eps
seps = seps(obi);
logp = -inf;
cov = zeros(nent, nent);
cov2 = zeros(nent, nent); %this covariance matrix is for gradient calculation

itr = 0;
diff = [];
psig = 0;
pseps = 0;
ph = 0;


cycle = 0;
while cycle < 1
    step_size = 0.001;
    itr = 0;
    %% training hyperparameter
    while itr<20
        %calculate K for this round of update
        for i = 1:nent
            for j = 1:nent
                [x1, y1] = unvec(i, nn);
                [x2, y2] = unvec(j, nn);
                dist = (x1-x2)^2+(y2-y1)^2;
                cov(i,j) = sigma^2*exp(-l*dist/2);
            end
        end
        origcov = cov;
        K = cov(obi, obi) + diag(seps.^2);
        L = chol(nearestSPD(K))';
        beta = L\filledValue(1:filled)';
        alpha = L'\beta;
        aat = alpha*alpha';
        invK = K\eye(filled);
        ami = aat - invK;

        %save old derivative
        oldpsig = psig;
        oldpseps = pseps;

        %calculate partial derivative for sigma
        psig = trace(ami*origcov(obi, obi))/sigma;

        %partial derivative for \sigma_\eps
        pseps = seps.*diag(ami);

        %update the hyperparameters 
        sigma = sigma + step_size*psig;
%         seps = seps + step_size*pseps;

        %calculate the log probability
        oldlogp = logp;
        % do not directly calculate determinant, it will blow up to inf
        % instead use cholesky decomposition
        % should have term - nent*log(2*pi)/2, but omit since it is a constant
        logp = -0.5*(beta'*beta)-0.5*sum(log(diag(L)));    
        fprintf('iter1 num %i, loss: %1.5e\n',itr,-logp);

        %implement backtracking line search to determine stepsize
        if -logp > -oldlogp - 0.5*step_size*norm([oldpsig; oldpseps])^2;
            step_size = step_size*0.5;
        end

        itr = itr+1;
        if logp == inf ||logp == -inf
            break
        end
    end

    %% train l separately
    itr = 0;
    step_size = 1;
    while itr<20
        %calculate K for this round of update
        for i = 1:nent
            for j = 1:nent
                [x1, y1] = unvec(i, nn);
                [x2, y2] = unvec(j, nn);
                dist = (x1-x2)^2+(y2-y1)^2;
                cov(i,j) = sigma^2*exp(-l*dist/2);
                cov2(i,j) = -sigma^2*dist*exp(-l*dist/2)/2;
            end
        end
        origcov = cov;
        K = cov(obi, obi) + diag(seps.^2);
        L = chol(nearestSPD(K))';
        beta = L\filledValue(1:filled)';
        alpha = L'\beta;
        aat = alpha*alpha';
        invK = K\eye(filled);
        ami = aat - invK;

        %save old derivative
        oldph = ph;

        %partial derivative for l
        pl = 0.5*trace(ami*cov2(obi,obi));
        ph = -2*pl/h^3;

        %update the hyperparameters
        h = h + step_size*ph;
        l = 1/h^2;

        %calculate the log probability
        oldlogp = logp;
        % do not directly calculate determinant, it will blow up to inf
        % instead use cholesky decomposition
        % should have term - nent*log(2*pi)/2, but omit since it is a constant
        logp = -0.5*(beta'*beta)-0.5*sum(log(diag(L)));
        fprintf('iter2 num %i, loss: %1.5e\n',itr,-logp);

        %implement backtracking line search to determine stepsize
        if -logp > -oldlogp 
            step_size = step_size*0.5;
        end

        itr = itr+1;
        if logp == inf ||logp == -inf
            break
        end
    end
cycle = cycle+1;
end

%% learning variance
%fill in learned variance
variance = zeros(dim);
for i = 1:filled
    [x,y] = unvec(obi(i), nn);
    variance(x,y) = seps(i);
end
% variance = DRM(variance,.5,.1,.1,.5,nn,5,2,1);
variance = vanilla_gpr(variance);

%% update covariance matrix
% note that we enumerate the entries by row, not column, caution when
% reshape
% cov = cov+diag(reshape(variance', nent,1).^2);
cov = cov+diag(reshape(dWvar(1:2:100, 1:2:100), nent,1));
