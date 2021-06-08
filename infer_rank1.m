%synthetic data 
nb = np*nn; %number of observations
% tdW = sChange; %true dW

%% Network parameters
k = 2;
Rate_Novel = gamrnd(k,meanRate/k,nn,np);
Rate_Novel(Rate_Novel>50) = 50;

BinConn = (rand(nn)<cd);    % Structural connectivity (it is assumed to be known)
if paramno==0
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
% h = h*10;
% h = h - mean(h)*mean(StrengthConn, 'all')*nn*cd/meanRate; % since \sum W\Del r is normally distributed


% using Gaussian Regression
% hyperparameter inference
model_selection_lbgfs

% Construct the giant covariance matrix, h is affine observation
%  hh   hx 
%  xh   xx

%upper right
hx = hxgen(sigma, l, vb_index, C, rf, nent); 

%upper left
K = Kgen(sigma, l, seps_neuron, vb_index, C, rf);
hh = K+diag(seps^2*ones(length(K),1));

%lower left
xh = hx';

K = hh;
Ks = hx;
% Kss = xx;


% calculate posterior
try
    L = chol(K)';
catch
    L = chol(nearestSPD(K))';
end
alpha = generalizedMinimalResidualMethod(K, h);
%     v = L\Ks;
%     K_pos = Kss - v'*v; %posterior variance 
mu_pos =Ks'*alpha;
ret = reshape(mu_pos, nf, nf);

function hx = hxgen(sigma, l, vb_index, C, rf, nent)
   nvb = length(vb_index);
   hx = zeros(nvb, nent);
   for i = 1:nvb
       for j = 1:nent
         npost1 = vb_index(i); % index of post-synaptic neuron 1
         rpost1 = rf(npost1); % firing rate of post-synaptic neuron 1
         rpre1 = rf(C(npost1,:));
         [rpre2, rpost2] = unvec(j, 50);
         for m = 1:length(rpre1)
             dist = (rpost1 - rpost2)^2+(rpre1(m) - rpre2)^2;
             hx(i,j) = hx(i,j) + sigma^2*exp(-dist/(2*l^2))*rpre1(m);
         end
       end
   end
end

function ret = isc(r1, r2, cut) %if in the same cluster
    ret = 0;
    for i = 1:length(cut)-1
        if(r1>=cut(i)&&r2>=cut(i)&&r1<cut(i+1)&&r2<cut(i+1))
            ret = 1;
        end
    end
end


