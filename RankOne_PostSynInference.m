clear all;clc;close all
%% Network parameters
Nneuron = 1000; % total number of neurons in the network
Npattern = 1;   % number of patterns it learned
meanRate = 10;  % Average firing rate before learning
k = 2;
Rate_Novel = gamrnd(k,meanRate/k,Nneuron,Npattern);

p = 0.1;        % probability of connection
BinConn = (rand(Nneuron)<p);    % Structural connectivity (it is assumed to be known)
StrengthConn = 1/(Nneuron*p)*(2*rand(Nneuron)+4);   % Strength of connection before learning
WRec_Novel = 1/meanRate*BinConn.*StrengthConn;

IExt = Rate_Novel-WRec_Novel*Rate_Novel;            % External input (assumed to be unchanged with learning)
%% Synaptic plasticity and firing rate after learning for each pattern
DelW_Strength = zeros(Nneuron,Nneuron,Npattern);
DelW_Strength_total = zeros(Nneuron);
for i = 1:Npattern
    DelW_Strength(:,:,i) = BinConn.*(1/(Nneuron*p)*.1/meanRate^3*(Rate_Novel(:,i).^2-meanRate*Rate_Novel(:,i))*(Rate_Novel(:,i)'-meanRate));
    DelW_Strength_total = DelW_Strength_total+DelW_Strength(:,:,i);
end

DelW = DelW_Strength_total;

WRec_Fam = WRec_Novel + DelW;

Rate_Fam = (eye(Nneuron) - WRec_Fam)\IExt;
%%  Inference on the post-synaptic dependence
Diff_Rate = Rate_Fam-Rate_Novel;
DelI_Rate = WRec_Novel*Diff_Rate;
% DelI_Rate_appr = 1/p*BinConn.*repmat(mean(WRec_Novel,2),1,Nneuron)*Diff_Rate;

Rate_Novel_vector = reshape(Rate_Novel,[],1);
Diff_Rate_vector = reshape(Diff_Rate,[],1);
[R,I] = sort(Rate_Novel_vector);
% post-synaptic dependence
figure;plot(R,Diff_Rate_vector(I),'.')

% no correlation between W*DelR and DelI (or DelR here)
ipattern = 1;
figure; plot(Diff_Rate(:,ipattern),DelI_Rate(:,ipattern),'.')
