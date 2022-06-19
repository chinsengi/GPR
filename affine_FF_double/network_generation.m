nn = 500;
paramno = 3;
cd = .1; % connection density
mxp = 1; % max # of patterns
nf = 50; % number of frequency levels
meanRate = 20;  % Average firing rate before learning
dim = [nf,nf];
nent = prod(dim);
sf = dim(1)/50; %scale factor
np = 1;
nb = np*nn; %number of observations

load(['paramset' num2str(paramno)]);
X = dW1(1:2:end, 1:2:end);
sChange = X;

rng(seed);
for i = 1:nf
    for j = 1:nf
        sChange(i,j) = sqrt(dWvar(i*(100/nf),j*(100/nf)))*randn()+X(i,j);
    end
end
sChange = (sChange - 0.5)/.5;

%% Network parameters
rng(seed);
k = 2;
Rate_Novel = gamrnd(k,meanRate/k,nn,np);
Rate_Novel(Rate_Novel>50) = 50;

BinConn = (rand(nn)<cd);    % Structural connectivity (it is assumed to be known)
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