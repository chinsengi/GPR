%% synthetic data 
% Network parameters
meanRate = 15;
k = 2;
npost = 100;
Rate_Pre = gamrnd(k,meanRate/k,nn,np);
% Rate_Pre = rand(nn, np)*50;
% Rate_Pre = lognrnd(2,0.5,nn,np);
Rate_Pre(Rate_Pre>70) = 70;
active_protocal = false;

BinConn = (rand(npost,nn)<cd);    % Structural connectivity (it is assumed to be known)
% StrengthConn = 1/(nn*cd)*(20*rand(npost,nn)+5);   % Strength of connection before learning
StrengthConn = lognrnd(-2.5,0.3,npost,nn);   % Strength of connection before learning
if expno>=9 && expno<=12
    StrengthConn = StrengthConn*5/(expno-8);
end


WFF_Novel = BinConn.*StrengthConn;
Rate_Novel = WFF_Novel*Rate_Pre;
c = min(Rate_Novel)-0.1;
Rate_Novel = Rate_Novel-c;
if expno>=9 && expno<= 12
    Rate_Novel = linspace(1, 50, length(Rate_Novel));
end
Rate_Novel(Rate_Novel>70) = 70;

if active_protocal
    high_thres = 50; % neurons firing at rate above this value are considered high firing rate neurons. 
    tmp = Rate_Pre(Rate_Pre>high_thres);
    ninh = ceil(numel(tmp));
    tmp(1:ninh) = rand(ninh, 1)*5+0.1;
    Rate_Pre(Rate_Pre>high_thres) = tmp;
    tmp = Rate_Novel(Rate_Novel>high_thres);
    nhigh = numel(tmp);
    tmp(1:nhigh) = rand(nhigh, 1)*40+0.1;
    Rate_Novel(Rate_Novel>high_thres) = tmp;
end

Rate_Pre = ceil(Rate_Pre);
Rate_Novel = ceil(Rate_Novel);
DelW = zeros(npost,nn);
tmp = [];
for i = 1:npost
    for j = 1:nn
        if BinConn(i,j)~=0
            tmp = [tmp; ceil(Rate_Novel(i)), Rate_Pre(j)];
            DelW(i,j) = sChange(Rate_Novel(i), Rate_Pre(j));
        end
    end
end 
WFF_Fam = WFF_Novel + DelW;

Rate_Fam = WFF_Fam*Rate_Pre-c;

rfpre = Rate_Pre;
rf = Rate_Novel;
C = BinConn;

% set missing connections
if expno==2
    tmp = C(BinConn~=0);
    C(BinConn~=0) = rand(size(tmp))< 0.9;
elseif expno==7
    tmp = C(BinConn~=0);
    C(BinConn~=0) = rand(size(tmp))< 0.8;
end

% add input noise
if expno>2 && expno < 7
    if expno==3 || expno == 4
        input_noise_level = 0.1;
    elseif expno==5 || expno==6
        input_noise_level = 0.5;
    end
    
    rf = rf + input_noise_level*randn(size(rf));
    rfpre = rfpre + input_noise_level*randn(size(rfpre));
end

% sample part of the observations
% vb stands for valid observation, nvb is the number of valid observations
vb_index = randperm(npost, nvb);
h = DelW*Rate_Pre;
h = h(vb_index,:);


% using Gaussian Regression
[param_1, bestlogp, postgrad] = train_param(model_selection_method, h, vb_index,...
                                                C, rf, rfpre, nent, input_noise);
exp(param_1)
if input_noise
    dW = infer(param_1, h, C, rf, rfpre, nent, vb_index, postgrad);
else
    dW = infer(param_1, h, C, rf, rfpre, nent, vb_index);
end
ret = reshape(dW, nf, nf);   

% calculate batch average posterior mean
function dW = infer(param, h, C, rf, rfpre, nent, vb_index, postgrad)
    [sigma, l, sn] = extract_param(param);
    if nargin>7
        hh = Kgen(param, vb_index, C, rf, rfpre, postgrad);
    else
        hh = Kgen(param, vb_index, C, rf, rfpre);
    end
    K = hh;
    Ks = hxgen(param, C, rf, rfpre, nent, vb_index); 
    alpha = generalizedMinimalResidualMethod(K, h);
    mu_pos =Ks'*alpha;
    dW = mu_pos;
end

function [param_true, bestlogp, postgrad] = train_param(model_selection_method, h, vb_index,...
                                                C, rf, rfpre, nent, input_noise)
    % initialize parameter
    sn_cand = [0.005];
    sig_cand = [0.02];
    postgrad = 0;
    
    bestlogp = -inf;
    for sn_ind = 1:length(sn_cand)
        for sig_ind = 1:length(sig_cand)
            sigma = sig_cand(sig_ind);
            l = 40;
            seps_neuron = sn_cand(sn_ind); % seps_neuron is the noise param for the original map (instead of the affine one)

            %current initial starting point. 
            X0 = log([sigma; l; seps_neuron]);
%             param_true = X0; return;
            % Construct the giant covariance matrix, h is affine observation
            % x is synaptic plasticity rule
            %  hh   hx 
            %  xh   xx
            % if there is input noise we want to obtain analytical gradient first 
            % so we need the following covariance matrix
            %  hh   hg
            %  gh   gg

            if input_noise
                hh = Kgen(X0, vb_index, C, rf, rfpre);
                K = hh;
                alpha = generalizedMinimalResidualMethod(K, h);
                % cov(h, grad)
                hg = hggen(X0, vb_index, C, rf, rfpre, nent*2);
                postgrad =hg'*alpha;
                postgrad = reshape(postgrad, 2, nent);
                %initial input noise level
                innoise = 0.1;
                X0 = [X0; log(innoise)];
            end

            % hyperparameter inference
            if input_noise
                [x_min, logp] = model_selection(X0, model_selection_method, h, vb_index,...
                                            C, rf, rfpre, postgrad);
            else
                [x_min, logp] = model_selection(X0, model_selection_method, h, vb_index,...
                                            C, rf, rfpre);
            end

            if logp > bestlogp
                param_true = x_min;
                bestlogp = logp;
            end
        end
    end
end

function hx = hxgen(param, C, rf, rfpre, nent, vb_index)
   [sigma, l] = extract_param(param);
   nvb = length(vb_index);
   hx = zeros(nvb, nent);
   for i = 1:nvb
       for j = 1:nent
         npost1 = vb_index(i); % index of post-synaptic neuron 1
         rpost1 = rf(npost1); % firing rate of post-synaptic neuron 1
         rfpre1 = rfpre(C(npost1,:));
         [rpre2, rpost2] = unvec(j, round(sqrt(nent)));
         for m = 1:length(rfpre1)
             dist = (rpost1 - rpost2)^2+(rfpre1(m) - rpre2)^2;
             hx(i,j) = hx(i,j) + sigma^2*exp(-dist/(2*l^2))*rfpre1(m);
         end
       end
   end
end

function hg = hggen(param, vb_index, C, rf, rfpre, nent)
    [sigma, l] = extract_param(param);
    nvb = length(vb_index);
    hg = zeros(nvb, nent);
    for i = 1:nvb
        for j = 1:2:nent
            npost1 = vb_index(i); % index of post-synaptic neuron 1
            rpost1 = rf(npost1); % firing rate of post-synaptic neuron 1
            rfpre1 = rfpre(C(npost1,:));
            [rpre2, rpost2] = unvec((j+1)/2, 50);
            for m = 1:length(rfpre1)
                dist = (rpost1 - rpost2)^2+(rfpre1(m) - rpre2)^2;
                hg(i,j) = hg(i,j) + sigma^2*exp(-dist/(2*l^2))*rfpre1(m)...
                                                     *(rfpre1(m)-rpre2)/l^2;
                hg(i,j+1) = hg(i,j+1)+sigma^2*exp(-dist/(2*l^2))*rfpre1(m)...
                                                     *(rpost1-rpost2)/l^2;
            end
        end
    end
end