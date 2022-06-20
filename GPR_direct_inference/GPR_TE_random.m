function [result,result2,result3,result4,rs,rs2,rs3,rs4,snr,snr2,M1,M2,M3,M4,M_sample,M_sample2]=GPR_TE_random(X,dWvar,heterosig,heterol,heteroseps,nsim,ntrial,srs)

%% driver file - everything start from here
f2=1:50;
[x2d,y2d]=meshgrid(f2,f2);
dsr=10;
nsr=length(srs);
nf = 50; %number of frequency quantum
dim = [nf,nf];
nent = prod(dim);
sf = dim(1)/nf; %scale factor
tdW = zeros(dim);
rho0=0.5;

noise=sqrt(dWvar);
dW=(X-rho0)/rho0;

%% active sampling and inference
result=zeros(ntrial,nsr,nsim);
result2=zeros(ntrial,nsr,nsim);
result3=zeros(ntrial,nsr,nsim);
result4=zeros(ntrial,nsr,nsim);

rs=zeros(ntrial,nsr,nsim);
rs2=zeros(ntrial,nsr,nsim);
rs3=zeros(ntrial,nsr,nsim);
rs4=zeros(ntrial,nsr,nsim);

snr=zeros(ntrial,nsr,nsim);
snr2=zeros(ntrial,nsr,nsim);
M1=cell(ntrial,nsr,nsim);
M2=cell(ntrial,nsr,nsim);
M3=cell(ntrial,nsr,nsim);
M4=cell(ntrial,nsr,nsim);
rng('default')
seedForSimulation = randi(10000,1,nsim);
M_sample=cell(ntrial,nsr,nsim);
M_sample2=cell(ntrial,nsr,nsim);
parfor kk=1:nsim
    rng(seedForSimulation(kk));
    [x2d,y2d]=meshgrid(f2,f2); 
    ndW=cell(ntrial,nsr,nsim);
    
    for jj=1:nsr
            [M0,indexinit]=voronoisampling(dim,srs(jj),X);
            indexKnown2=indexinit;
            indexKnown=randperm(2500,srs(jj))';
            ndW_all=zeros(size(dW));
            for ii=1:ntrial

                tdW=reshape(X(:)+noise(:).*randn(size(X(:))),size(X));
                tdW = (tdW - rho0)/rho0;
                ndW{ii,jj,kk}=tdW;
                ndW_all=ndW_all+ndW{ii,jj,kk};
                ndW_mean=ndW_all/ii;
                observed=zeros(size(dW));
                observed(indexKnown)=ndW_mean(indexKnown);
                observed2=zeros(size(dW));
                observed2(indexKnown2)=ndW_mean(indexKnown2);
%                 pred3=svd_complete(observed,observed);
%                 result3(ii,jj,kk)=norm(pred3(:)- dW(:))/norm(dW(:));
%                 rs3(ii,jj,kk)=sum((pred3(:)-dW(:)).^2)/sum((dW(:)-mean(dW(:))).^2);%rsquare(X(:),pred3(:));
%            
                M_sample{ii,jj,kk}=observed;
                M_sample2{ii,jj,kk}=observed2; %sample for voronoi
                
                [result(ii,jj,kk),rs(ii,jj,kk),reconstW]=GPR_core(observed,dW,indexKnown,length(indexKnown),observed(indexKnown)',heterosig,heterol,heteroseps);    
                M1{ii,jj,kk}=reconstW(:);
                
                xobserved=x2d(indexKnown);
                yobserved=y2d(indexKnown);
                zobserved=ndW_mean(indexKnown);
                fitresult=TaylorExpansionFit(xobserved,yobserved,zobserved);
                pred2=fitresult(x2d,y2d);
                %pred2=SynapticLinearRegression(indexKnown',observed(indexKnown),2,4);
                M2{ii,jj,kk}=pred2(:);
                result2(ii,jj,kk)=norm(pred2(:)- dW(:))/norm(dW(:));
                rs2(ii,jj,kk)=sum((pred2(:)-dW(:)).^2)/sum((dW(:)-mean(dW(:))).^2);
                
                
                %%%%%%%%%%GPR
                [result3(ii,jj,kk),rs3(ii,jj,kk),reconstW3]=GPR_core(observed2,dW,indexKnown2,length(indexKnown2),observed2(indexKnown2)',heterosig,heterol,heteroseps);    
                M3{ii,jj,kk}=reconstW3(:);                
                
                xobserved=x2d(indexKnown2);
                yobserved=y2d(indexKnown2);
                zobserved=ndW_mean(indexKnown2);
                fitresult2=TaylorExpansionFit(xobserved,yobserved,zobserved);
                pred4=fitresult2(x2d,y2d);
                %pred2=SynapticLinearRegression(indexKnown',observed(indexKnown),2,4);
                M4{ii,jj,kk}=pred4(:);
                result4(ii,jj,kk)=norm(pred4(:)- dW(:))/norm(dW(:));
                rs4(ii,jj,kk)=sum((pred4(:)-dW(:)).^2)/sum((dW(:)-mean(dW(:))).^2);

                                
                istnoise=ndW_mean-dW;
                snr(ii,jj,kk)=norm(istnoise(:))/norm(dW(:));
                %snr(ii,jj,kk)=mean(istnoise(indexKnown));
                snr2(ii,jj,kk)=norm(istnoise(indexKnown))/norm(dW(indexKnown));
            end
    end
end

end
