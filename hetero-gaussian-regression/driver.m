%% driver file - everything start from here
nf = 50; %number of frequency quantum
dim = [nf,nf];
nent = prod(dim);
sf = dim(1)/nf; %scale factor
tdW = zeros(dim);
initparam = true;

rule = 0; % 0 stands for calcium-based rule, 1 stands for BCM based rule
if rule == 1
    theta1 = 0;
    theta2 = 8; %BCM rule
    for i = 1:nf
        g = i/sf - 25;
        for j = 1:nf
            if i <= theta1
                f = 0;
            else
                f = (j/sf-theta1)*(j/sf-theta2);
            end
            tdW(i,j) = f*g;
        end
    end
    tdW = tdW/10000;
    X = tdW;
    folder = "figs_bcm";
elseif rule == 0
    load(['paramset' num2str(paramno)]);
    X = dW1(2:2:end, 2:2:end);
    dWvar = dWvar(1:2:end, 1:2:end);
%     X = dWs{paramno};
%     dWvar = dWvars{paramno};
    tdW = X;
    tdW(isnan(tdW)) = 1;
    dim = size(tdW);
    zfolder = "figs_cal";
else
    load('prMatrix');
    tdW = X-1e-6;
    dim = size(X);
end

%% add noise to data
if noise_index == 0
    % add uniform noise
    tdW = X+randn(size(X));
else 
    % add structured noise
%     seed = randi(10000);
%     save_seed(paramno+1, trial) = seed;
%     save('./result/save_seed', 'save_seed');
%     seed = 78819;
%     seed = save_seed(paramno+1, trial);
%     rng(seed);
    for i = 1:nf
        for j = 1:nf
            tdW(i,j) = sqrt(dWvar(i,j))*randn()+tdW(i,j);
        end
    end
end

% normalization with rho0 = 0.5
rho0 = 0.5;
tdW = (tdW - rho0)/rho0;
X = (X-rho0)/rho0;

%% generate latin grid
% init = randperm(nn);
% init = flip(1:dim(1)); %para-diagonal

Gaussian_regression_2
