%% synthetic data 
% Network parameters
nn = 500;
if expno ==1 
    nnE = 400;
    nnI = nn - nnE;
elseif expno > 1
    nnE = round(nn/2);
    nnI = nn - nnE;
end
npost = 100;
k = 2;
meanRateE = 15;
RateE_Pre = gamrnd(k,meanRateE/k,nnE,np);
RateE_Pre(RateE_Pre>70) = 70;
RateE_Pre = ceil(RateE_Pre);

meanRateI = 20; % 20
RateI_Pre = gamrnd(k,meanRateI/k,nnI,np);
RateI_Pre(RateI_Pre>70) = 70;
RateI_Pre = ceil(RateI_Pre); 

BinConn = (rand(npost,nn)<cd);    % Structural connectivity (it is assumed to be known)
% StrengthConn = 1/(nn*cd)*(20*rand(npost,nn)+5);   % Strength of connection before learning
StrengthConnE = lognrnd(-2.5,0.3,npost,nnE);   % Strength of connection before learning
StrengthConnI = lognrnd(-2.5,0.3,npost,nnI);   % Strength of connection before learning

WEFF_Novel = BinConn(:,1:nnE).*StrengthConnE;
WIFF_Novel = BinConn(:,nnE+1:nn).*StrengthConnI;

Rate_Novel = WEFF_Novel*RateE_Pre - WIFF_Novel*RateI_Pre;
c = min(Rate_Novel) - 0.1;
Rate_Novel = Rate_Novel - c;
Rate_Novel(Rate_Novel<1) = 1;
Rate_Novel(Rate_Novel>70) = 70;
if expno>1
    Rate_Novel = linspace(1, 50, length(Rate_Novel));
end
Rate_Novel = ceil(Rate_Novel);

DelWE = zeros(npost,nnE);
pdE = [];
for i = 1:npost
    for j = 1:nnE
        if BinConn(i,j)~=0
            pdE = [pdE; Rate_Novel(i), RateE_Pre(j)];
            DelWE(i,j) = sChange(Rate_Novel(i), RateE_Pre(j));
        end
    end
end 
WEFF_Fam = WEFF_Novel + DelWE;
WEFF_Fam(WEFF_Fam<0) = 0;

R0 = meanRateI;
DelWI = zeros(npost,nnI);
tmp = 1:nf; 
if expno == 1
    scaleFactor = 80000; % set the scaling factor to be 80000 to have same magnitude
elseif expno == 2
    scaleFactor = 160000;
elseif expno == 3
    scaleFactor = 40000;
elseif expno == 4
    scaleFactor = 80000*2/3;
elseif expno == 5
    scaleFactor = 120000;
elseif expno == 6
    scaleFactor = -80000;
elseif expno == 7
    scaleFactor = 80000*4;
elseif expno == 8
    scaleFactor = 80000*3;
elseif expno == 9
    scaleFactor = 80000;
elseif expno == 10 
    scaleFactor = 80000/3;
elseif expno == 11
    scaleFactor = 80000/4;
end
XI = (tmp'-R0)*tmp/scaleFactor;
stdInoise = std(sChange(:) - X(:))*80000/scaleFactor; 
pdI = [];
for i = 1:npost
    for j = 1:nnI
        if BinConn(i,j+nnE)~=0
            pdI = [pdI; Rate_Novel(i) RateI_Pre(j)];
            DelWI(i,j) = (Rate_Novel(i)-R0)*RateI_Pre(j)/scaleFactor ...
                                                         +randn()*stdInoise;
%             DelWI(i,j) = sChange(Rate_Novel(i), RateI_Pre(j));
        end
    end
end 
WIFF_Fam = WIFF_Novel + DelWI;
WIFF_Fam(WIFF_Fam<0) = 0;
Rate_Fam = WEFF_Fam*RateE_Pre - WIFF_Fam*RateI_Pre - c;
Rate_Fam(Rate_Fam<0) = 0;

rfpreI = RateI_Pre;
rfpreE = RateE_Pre;
rf = Rate_Novel;
C = BinConn;

% sample part of the observations
% vb stands for valid observation, nvb is the number of valid observations
vb_index = randperm(npost, nvb);
% vb_index = sparse_sampling(nvb, nb, C);
% h = DelW*Rate_Fam + randn(nb,1)*.6647 +  4.1530;
% h = DelW*Rate_Fam + WRec_Novel*Diff_Rate;
% assert(norm(Diff_Rate-h)<1e-10);
h = DelWE*RateE_Pre - DelWI*RateI_Pre;
h = h(vb_index,:);

%calculate the inference range
% pairsE = []; pairsI = [];
% for i = 1:length(vb_index)
%     pairsE = [pairsE; pdE{vb_index(i)}];
%     pairsI = [pairsI; pdI{vb_index(i)}];
% end
% center = mean(pairsE,1);
% stdx = std(pairsE(:,1));
% stdy = std(pairsE(:,2));
% irEx = max(1,ceil(center-stdx)):min(50, ceil(center+stdx)); 
% irEy = max(1,ceil(center-stdy)):min(50, ceil(center+stdy));
% 
% center = mean(pairsI,1);
% stdx = std(pairsI(:,1));
% stdy = std(pairsI(:,2));
% irIx = max(1, ceil(center(1)-stdx)):min(50, ceil(center+stdx)); 
% irIy = max(1, ceil(center(1)-stdy)):min(50, ceil(center+stdy));



% using Gaussian Regression
[param_1, bestlogp] = train_param(model_selection_method, h, vb_index,...
                            C, rf, rfpreI, rfpreE, input_noise);
exp(param_1)
if input_noise
    dW = infer(param_1, h, C, rf, rfpreI, rfpreE, nent, vb_index,...
                                                      postgrad);
else
    dW = infer(param_1, h, C, rf, rfpreI, rfpreE, nent, vb_index);
end
dWI = reshape(dW(1:nent), nf, nf);   
dWE = reshape(dW(nent+1:end), nf, nf);

% calculate batch average posterior mean
function dW = infer(param, h, C, rf, rfpreI, rfpreE, nent, vb_index,...
                                                                  postgrad)
    nvb = numel(vb_index);
    [sigmaE, lE, snE, sigmaI, lI, snI, innoise] = extract_param(param);
    if nargin>8
        hh = Kgen(param, vb_index, C, rf, rfpreI, rfpreE, postgrad);
    else
        hh = Kgen(param, vb_index, C, rf, rfpreI, rfpreE);
    end
    K = hh;
    Ks = hxgen(param, C, rf, rfpreI, rfpreE, nent, vb_index); 
    alpha = generalizedMinimalResidualMethod(K, h);
    mu_pos =Ks'*alpha;
    dW = mu_pos;
end

function [param_true, bestlogp] = train_param(model_selection_method, h, vb_index,...
                                   C, rf, rfpreI, rfpreE, input_noise)
    % initialize parameter
    sn_cand = [0.002];
    sig_cand = [0.01];

    bestlogp = -inf;
    for sn_ind = 1:length(sn_cand)
        for sig_ind = 1:length(sig_cand)
            sigma = sig_cand(sig_ind);
            l = 40;
            seps_neuron = sn_cand(sn_ind); % seps_neuron is the noise param for the original map (instead of the affine one)

            %current initial starting point. 
            X0 = log([sigma; l; seps_neuron]);
            X0 = [X0;X0];

            % Construct the giant covariance matrix, h is affine observation
            % x is synaptic plasticity rule
            %  hh   hx 
            %  xh   xx
            % if there is input noise we want to obtain analytical gradient first 
            % so we need the following covariance matrix
            %  hh   hg
            %  gh   gg

            if input_noise
                hh = Kgen(X0, vb_index, C, rf, rfpreI, rfpreE);
                hh = hh+diag(seps^2*ones(length(hh),1));
                K = hh;
                alpha = generalizedMinimalResidualMethod(K, h-noise_mu);
                % cov(h, grad)
                hg = hggen(sigma, l, vb_index, C, rf, rfpreI, rfpreE, nent*2);
                postgrad =hg'*alpha;
                % first row of postgrad is derivative w.r.t. x,
                % second row of postgrad is derivative w.r.t. y;
                postgrad = reshape(postgrad, 2, nent); 
                %initial input noise level
                innoise = 0.1;
                X0 = [X0 log(innoise)];
            end

            % hyperparameter inference
            if input_noise
                [x_min, logp] = model_selection(X0,...
                  model_selection_method, h, vb_index,C, rf, ...
                  rfpreI, rfpreE, postgrad);
            else
                [x_min, logp] = model_selection(X0,...
                    model_selection_method, h, vb_index,...
                              C, rf, rfpreI, rfpreE);
            end

            if logp > bestlogp
                param_true = x_min;
                bestlogp = logp;
            end
        end
    end
end

function hx = hxgen(param, C, rf, rfpreI, rfpreE, nent, vb_index)
   nvb = numel(vb_index);
   nnE = numel(rfpreE);
   nnI = numel(rfpreI);
   nn = nnE+nnI;
   [sigmaE, lE, snE, sigmaI, lI, snI, innoise] = extract_param(param);
   hx = zeros(nvb, nent*2);
   for i = 1:nvb
       for j = 1:nent
         npost1 = vb_index(i); % index of post-synaptic neuron 1
         rpost1 = rf(npost1); % firing rate of post-synaptic neuron 1
         rpre1I = rfpreI(C(npost1,nnE+1:nn));
         [rpre2, rpost2] = unvec(j, 50);
         for m = 1:length(rpre1I)
             dist = (rpost1 - rpost2)^2+(rpre1I(m) - rpre2)^2;
             hx(i,j) = hx(i,j) - sigmaI^2*exp(-dist/(2*lI^2))*rpre1I(m);
         end
       end
   end
   for i = 1:nvb
       for j = 1:nent
         npost1 = vb_index(i); % index of post-synaptic neuron 1
         rpost1 = rf(npost1); % firing rate of post-synaptic neuron 1
         rpre1E = rfpreE(C(npost1, 1:nnE));
         [rpre2, rpost2] = unvec(j, 50);
         for m = 1:length(rpre1E)
             dist = (rpost1 - rpost2)^2+(rpre1E(m) - rpre2)^2;
             hx(i,j+nent) = hx(i,j+nent) + sigmaE^2*exp(-dist/(2*lE^2))*rpre1E(m);
         end
       end
   end
end

function hg = hggen(sigma, l, vb_index, C, rf, rfpreI, rfpreE, nent)
    nvb = length(vb_index);
    nnE = numel(rfpreE);
    nnI = numel(rfpreI);
    nn = nnE+nnI;
    hg = zeros(nvb, nent);
    for i = 1:nvb
        for j = 1:2:nent
            npost1 = vb_index(i); % index of post-synaptic neuron 1
            rpost1 = rf(npost1); % firing rate of post-synaptic neuron 1
            rpre1I = rfpreI(C(npost1,nnE+1:nn));
            rpre1E = rfpreE(C(npost1, 1:nnE));
            [rpre2, rpost2] = unvec((j+1)/2, 50);
            for m = 1:length(rpre1I)
                dist = (rpost1 - rpost2)^2+(rpre1I(m) - rpre2)^2;
                hg(i,j) = hg(i,j) + sigma^2*exp(-dist/(2*l^2))*rpre1I(m)...
                                                     *(rpre1(m)-rpre2)/l^2;
                hg(i,j+1) = hg(i,j+1)+sigma^2*exp(-dist/(2*l^2))*rpre1I(m)...
                                                     *(rpost1-rpost2)/l^2;
            end
            for m = 1:length(rpre1E)
                dist = (rpost1 - rpost2)^2+(rpre1E(m) - rpre2)^2;
                hg(i,j) = hg(i,j) + sigma^2*exp(-dist/(2*l^2))*rpre1E(m)...
                                                     *(rpre1(m)-rpre2)/l^2;
                hg(i,j+1) = hg(i,j+1)+sigma^2*exp(-dist/(2*l^2))*rpre1E(m)...
                                                     *(rpost1-rpost2)/l^2;
            end
        end
    end
end