%% initialize parameter
sigma = 1;
l = 0.0007;
hprime = sqrt(1/l);
seps = 0.2*ones(nent,1); %\sigma_\eps
logp = -inf;
cov = zeros(nent, nent);
cov2 = zeros(nent, nent); %this matrix is for gradient calculation

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
    while itr<40
        %calculate K for this round of update
%         for i = 1:nent
%             for j = 1:nent
%                 [x1, y1] = unvec(i, nn);
%                 [x2, y2] = unvec(j, nn);
%                 dist = (x1-x2)^2+(y2-y1)^2;
%                 cov(i,j) = sigma^2*exp(-l*dist/2);
%             end
%         end
        % magic mex
        cov = covgen(sigma, l, cov);
        origcov = cov;
        cov = cov + diag(seps.^2);
        K = A*cov*A';
        origK = A*origcov*A';
        L = chol(nearestSPD(K))';
        beta = L\h;
        alpha = L'\beta;
        aat = alpha*alpha';
        invL = inv(L);
        invK = invL'*invL;
        ami = aat - invK;

        %save old derivative
        oldpsig = psig;
        oldpseps = pseps;

        %calculate partial derivative for sigma
        psig = trace(ami*origK)/sigma;

        %partial derivative for \sigma_\eps
        pseps = zeros(nent, 1);
        for i = 1:nent
            tmp = A(:,i)*A(:,i)';
            pseps(i) = seps(i)*trace(ami*tmp);
        end
        
        %update the hyperparameters 
        sigma = sigma + step_size*psig;
        seps = seps + step_size*pseps;

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
    step_size = 10;
    while itr<40
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
        cov = cov + diag(seps.^2);
        K = A*cov*A';
        origK = A*cov2*A';
        L = chol(nearestSPD(K))';
        beta = L\h;
        alpha = L'\beta;
        aat = alpha*alpha';
        invL = inv(L);
        invK = invL'*invL;
        ami = aat - invK;

        %save old derivative
        oldph = ph;

        %partial derivative for l
        pl = 0.5*trace(ami*origK);
        ph = -2*pl/hprime^3;

        %update the hyperparameters
        hprime = hprime + step_size*ph;
        l = 1/hprime^2;

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
for i = 1:nent
    [x,y] = unvec(i, nf);
    variance(x,y) = seps(i);
end

%% update covariance matrix
% note that we enumerate the entries by row, not column, caution when
% reshape
cov = cov+diag(reshape(variance', nent,1).^2);
% cov = cov+diag(reshape(dWvar(1:2:100, 1:2:100), nent,1));
