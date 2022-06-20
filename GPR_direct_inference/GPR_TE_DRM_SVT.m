function [result,result2,result3,result4,rs,rs2,rs3,rs4,snr,snr2,M1,M2,M_sample]=GPR_TE_DRM_SVT(X,dWvar,heterosig,heterol,heteroseps,nsim,ntrial,srs)

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
opts.MAX_ITERS=10000;
opts.GEN_PLOTS=false;
opts.QUIET = true;
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
M_sample=cell(ntrial,nsr,nsim);
rng('default')
seedForSimulation = randi(10000,1,nsim);
for kk=1:nsim
    rng(seedForSimulation(kk));
    [x2d,y2d]=meshgrid(f2,f2); 
    ndW=cell(ntrial,nsr,nsim);
    
    for jj=1:nsr
 
            [M0,indexinit]=voronoisampling(dim,srs(jj),X);

            indexKnown=indexinit;
            %indexKnown=randperm(2500,srs(jj))';
            ndW_all=zeros(size(dW));
            for ii=1:ntrial
            %     stream = RandStream('dsfmt19937','Seed',1+i);
            %     ndW=reshape(dW(:)+noise(:).*randn(stream,size(dW(:))),size(dW));
                tdW=reshape(X(:)+noise(:).*randn(size(X(:))),size(X));
                tdW = (tdW - rho0)/rho0;
                ndW{ii,jj,kk}=tdW;
                ndW_all=ndW_all+ndW{ii,jj,kk};
                ndW_mean=ndW_all/ii;
                observed=zeros(size(dW));
                observed(indexKnown)=ndW_mean(indexKnown);
                
%                 pred3=svd_complete(observed,observed);
%                 result3(ii,jj,kk)=norm(pred3(:)- dW(:))/norm(dW(:));
%                 rs3(ii,jj,kk)=sum((pred3(:)-dW(:)).^2)/sum((dW(:)-mean(dW(:))).^2);%rsquare(X(:),pred3(:));
%                 
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
                
%                 result2(ii,jj,kk)=norm(pred2(indexKnown)- dW(indexKnown))/norm(dW(indexKnown));
%                 rs2(ii,jj,kk)=sum((pred2(indexKnown)-dW(indexKnown)).^2)/sum((dW(indexKnown)-mean(dW(indexKnown))).^2);
                pred3=DRM(observed*10000, 0.5, 0.5,0.1,0.06,50,8,2,1,opts)/10000;
                result3(ii,jj,kk)=norm(pred3(:)- dW(:))/norm(dW(:));
                rs3(ii,jj,kk)=sum((pred3(:)-dW(:)).^2)/sum((dW(:)-mean(dW(:))).^2);
                
                M_sample{ii,jj,kk}=observed;
                
                A=observed;
                B=A;
                B(B~=0)=1;
                lamnbda_tol = 0.01;
                tol = 1e-8;
                N = 1000;
                scale2=1;
                pred4 = MatrixCompletion(observed*scale2, B,N, 'nuclear', lamnbda_tol, tol, 0);
                pred4=pred4/scale2;
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
