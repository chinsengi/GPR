%synthetic data 
nb = np*nn; %number of observations

%% Network parameters
k = 2;
Rate_Novel = gamrnd(k,meanRate/k,nn,np);
Rate_Novel(Rate_Novel>50) = 50;

BinConn = (rand(nn)<cd);    % Structural connectivity (it is assumed to be known)
StrengthConn = 1/(nn*cd)*(2*rand(nn)+4);   % Strength of connection before learning
WRec_Novel = 1/meanRate*BinConn.*StrengthConn;
if paramno==0
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
    Rate_Fam(Rate_Fam<0) = 0;
    %%  Inference on the post-synaptic dependence
    Diff_Rate = Rate_Fam-Rate_Novel;
end


rf = Rate_Fam;
C = BinConn;

% add input noise
if input_noise
    input_noise_level = 0.1;
    rf = rf + input_noise_level*randn(size(rf));
end
% input_noise_level = 0.1;
% rf = rf + input_noise_level*randn(size(rf));

% sample part of the observations
% vb stands for valid observation, nvb is the number of valid observations
vb_index = randperm(nb);
vb_index = vb_index(1:nvb);
h = DelW*Rate_Fam + WRec_Novel*Diff_Rate;
h = h(vb_index,:);
% h = h*10;
% h = h - mean(h)*mean(StrengthConn, 'all')*nn*cd/meanRate; % since \sum W\Del r is normally distributed


% using Gaussian Regression
% initialize parameter
sn_cand = [0.1];
nm_cand = [1, 2];
seps_cand = [0.1, 0.2];
sig_cand = [0.02, 0.05];
bestlogp = -inf;
for sn_ind = 1:length(sn_cand)
    for nm_ind = 1:length(nm_cand)
        for seps_ind = 1:length(seps_cand)
            for sig_ind = 1:length(sig_cand)
                noise_mu = nm_cand(nm_ind);
                % sigma = sqrt(abs(mean(h) - noise_mu)/(nn*cd*meanRate));
                sigma = sig_cand(sig_ind);
                l = 40;
                seps = seps_cand(seps_ind);
                if heteroseps
                    seps = seps*ones(nvb,1); %\sigma_\epsilon
                end
                seps_neuron = sn_cand(sn_ind); % seps_neuron is the noise param for the original map (instead of the affine one)

                % Construct the giant covariance matrix, h is affine observation
                % x is synaptic plasticity rule
                %  hh   hx 
                %  xh   xx
                % if there is input noise we want to obtain analytical gradient first 
                % so we need the following covariance matrix
                %  hh   hg
                %  gh   gg

                if input_noise
                    hh = Kgen(sigma, l, seps_neuron, vb_index, C, rf);
                    hh = hh+diag(seps^2*ones(length(hh),1));
                    K = hh;
                    alpha = generalizedMinimalResidualMethod(K, h-noise_mu);
                    % cov(h, grad)
                    hg = hggen(sigma, l, vb_index, C, rf, nent*2);
                    postgrad =hg'*alpha;
                    postgrad = reshape(postgrad, 2, nent);
                    %initial input noise level
                    innoise = 0.1;
                end

                % hyperparameter inference
                model_selection_crossval

                if logp > bestlogp
                    param_true = x_min;
                    bestlogp = logp;
                end
            end
        end
    end
end
[sigma, l, seps, noise_mu, seps_neuron, innoise] = extract_param(param_true, nvb, heteroseps);

% calculate posterior mean
if input_noise
    hh = Kgen(sigma, l, seps_neuron, vb_index, C, rf, innoise, postgrad);
else
    hh = Kgen(sigma, l, seps_neuron, vb_index, C, rf);
end
if heteroseps
    hh = hh + diag(seps.^2);
else
    hh = hh+diag(seps^2*ones(length(hh),1));
end
K = hh;
Ks = hxgen(sigma, l, vb_index, C, rf, nent); 
alpha = generalizedMinimalResidualMethod(K, h-noise_mu);
% L = chol(nearestSPD(K))';
% alpha = L'\(L\h);
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

function hg = hggen(sigma, l, vb_index, C, rf, nent)
    nvb = length(vb_index);
    hg = zeros(nvb, nent);
    for i = 1:nvb
        for j = 1:2:nent
            npost1 = vb_index(i); % index of post-synaptic neuron 1
            rpost1 = rf(npost1); % firing rate of post-synaptic neuron 1
            rpre1 = rf(C(npost1,:));
            [rpre2, rpost2] = unvec((j+1)/2, 50);
            for m = 1:length(rpre1)
                dist = (rpost1 - rpost2)^2+(rpre1(m) - rpre2)^2;
                hg(i,j) = hg(i,j) + sigma^2*exp(-dist/(2*l^2))*rpre1(m)...
                                                     *(rpre1(m)-rpre2)/l^2;
                hg(i,j+1) = hg(i,j+1)+sigma^2*exp(-dist/(2*l^2))*rpre1(m)...
                                                     *(rpost1-rpost2)/l^2;
            end
        end
    end
end

function [sigma, l, seps, noise_mu, sn, innoise] = extract_param(X, nvb, heteroseps)
    % it's a hack, so that we don't need to know if there is input noise
    X = [X; 0];
    X = exp(X);
    if heteroseps
        sigma = X(1); l = X(2);
        seps = X(3:2+nvb);
        [noise_mu, sn, innoise] = feval(@(x) x{:}, num2cell(X(3+nvb:end)));
    else
        [sigma, l, seps, noise_mu, sn, innoise] = feval(@(x) x{:}, num2cell(X));
    end
end