function [result,rs,reconstW ]=GPR_core(dW,X,obi,filled,filledValue,heterosig,heterol,heteroseps)
dim=[50,50];
nent=2500;
nf=50;

dataMean = 0;
bestlogp = -inf;
sigma_guess = [0.1, 0.2, 0.05];
l_guess = [30 40 50];
seps_guess = [0.05 0.1 0.3];
for sigma_index = 1:3
    for l_index = 1:3
        for seps_index = 1:3
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
    K = K+diag(seps.^2);
else
    K = K+diag(seps(1)^2*ones(length(K),1));
end
alpha = generalizedMinimalResidualMethod(K, filledValue');
%         v = L\Ks;
%     K_pos = Kss - v'*v; %posterior variance
mu_pos = Ks'*alpha;
reconstW = dataMean + reshape(mu_pos, size(dW));

result= norm(reconstW(:)- X(:))/norm(X(:));
rs=sum((reconstW(:)-X(:)).^2)/sum((X(:)-mean(X(:))).^2);

% result= norm(reconstW(filled)- X(filled))/norm(X(filled));
% rs=sum((reconstW(filled)-X(filled)).^2)/sum((X(filled)-mean(X(filled))).^2);

result_logp = bestlogp;