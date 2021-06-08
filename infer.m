%synthetic data 

% rng(1);

nb = np*nn; %number of observations
tdW = sChange; %true dW
% tdW = ones(dim);
%     dW = zeros(dim); % dW with missing value
C = ((rand(nn,nn)/cd)<1); % random connection matrix
% try clustering
% if nn>=250
%     parse = randi([1 nn], 25,1);
%     parse = sort(parse);
%     parse = [1 parse' nn+1];
%     s = 0;
%     for i = 1:length(parse)-1
%         start = parse(i);
%         End = parse(i+1);
%         C(start:End-1, start:End-1) = ((rand(End-start,End-start)/cd)<5);
%         s=s+(End-start)^2;
%     end 
% end
% E = (s*cd*5+(nn*nn-s)*cd*.7)/(nn*nn);

%     rf = zeros(nn,np);
% for i = 1:np
%     mu = ceil(gamrnd(2,10, nn,1));
%     mu(mu>50) = 50;
%     sigma = rand(nn,nn)*2;
%     sigma = sigma*sigma';
%     sigma = sigma*3;
%     rf(:,i) = mvnrnd(mu, sigma);
% end

rf = gamrnd(2,10, nn,np); %synthetic firing rate for multiple pattern
rf(rf>50) = 50;
rf(rf<1) = 1;

%generate A matrix
A = sparse(nb, prod(dim));

for i = 1:np
    for j = 1:nn % process j-th neuron
        rowi = (i-1)*nn+j;
        ri = ceil(rf(j,i)*sf);
        rjs = C(j,:).*rf(:,i)';
        rjs = rjs(rjs~=0);
        for k = 1:length(rjs)
           coli = (ri-1)*dim(1)+ceil(rjs(k)*sf);
           A(rowi, coli) = A(rowi,coli) + rjs(k);
        end
    end
end

% sample part of the observations
% vb stands for valid observation, nvb is the number of it
vb_index = randperm(nb);
vb_index = vb_index(1:nvb);
A = A(vb_index, :);



h = A*tdW(:);

% using Gaussian Regression
% Construct the giant covariance matrix h is affine observation
%  hh   hx 
%  xh   xx

% lower right
% xx = zeros(nent, nent);
% for i = 1:nent
%     for j = i:nent
%         [x1, y1] = unvec(i, dim(1));
%         [x2, y2] = unvec(j, dim(1));
%         xx(i,j) = sigma^2*exp(-l*((x1-x2)^2+(y2-y1)^2)/2);
%     end
% end
% xx = xx+xx'- diag(diag(xx));
model_selection_structure2
xx = cov;

%upper right
hx = A*xx; 

%upper left
hh = hx*A';

%lower left
xh = hx';

K = hh;
Ks = hx;
Kss = xx;


% calculate posterior
% L = chol(nearestSPD(K))';
L = chol(nearestSPD(K))';
alpha = L'\(L\h);
%     v = L\Ks;
%     K_pos = Kss - v'*v; %posterior variance 
mu_pos =Ks'*alpha;
ret = reshape(mu_pos, nf, nf);

function ret = isc(r1, r2, cut) %if in the same cluster
    ret = 0;
    for i = 1:length(cut)-1
        if(r1>=cut(i)&&r2>=cut(i)&&r1<cut(i+1)&&r2<cut(i+1))
            ret = 1;
        end
    end
end

