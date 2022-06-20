function [result,result2,result3,rs,rs2,rs3,snr,snr2,M1,M2]=GPR_sampling(X,dWvar,heterosig,heterol,heteroseps,nsim,ntrial,srs)

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
rng('default')
seedForSimulation = randi(10000,1,nsim);
parfor kk=1:nsim
    rng(seedForSimulation(kk));
    [x2d,y2d]=meshgrid(f2,f2); 
    ndW=cell(ntrial,nsr,nsim);
    
    
    for jj=1:nsr
            [M0,indexinit]=voronoisampling(dim,srs(jj),X);
            indexKnown=indexinit;

            %[~,indexKnown2]=diagonalsampling(dim(1),srs(jj),X);

            indexKnown3=randperm(2500,srs(jj))';

            %indexKnown4=gridsampling(dim(1),srs(jj));

            %indexKnown5=boundsampling(dim(1),srs(jj))';
            
            ndW_all=zeros(size(dW));
            for ii=1:ntrial
            %     stream = RandStream('dsfmt19937','Seed',1+i);
            %     ndW=reshape(dW(:)+noise(:).*randn(stream,size(dW(:))),size(dW));
                tdW=reshape(X(:)+noise(:).*randn(size(X(:))),size(X));
                tdW = (tdW - rho0)/rho0;
                ndW{ii,jj,kk}=tdW;
                ndW_all=ndW_all+ndW{ii,jj,kk};
                ndW_mean=ndW_all/ii;
                %ndW_mean=tdW; %No average
                observed=zeros(size(dW));
                observed(indexKnown)=ndW_mean(indexKnown);
%                 observed2=zeros(size(dW));
%                 observed2(indexKnown2)=ndW_mean(indexKnown2);
                observed3=zeros(size(dW));
                observed3(indexKnown3)=ndW_mean(indexKnown3);
%                 observed4=zeros(size(dW));
%                 observed4(indexKnown4)=ndW_mean(indexKnown4);
%                 observed5=zeros(size(dW));
%                 observed5(indexKnown5)=ndW_mean(indexKnown5);
%                 pred3=svd_complete(observed,observed);
%                 result3(ii,jj,kk)=norm(pred3(:)- dW(:))/norm(dW(:));
%                 rs3(ii,jj,kk)=sum((pred3(:)-dW(:)).^2)/sum((dW(:)-mean(dW(:))).^2);%rsquare(X(:),pred3(:));               
                [result(ii,jj,kk),rs(ii,jj,kk),reconstW]=GPR_core(observed,dW,indexKnown,length(indexKnown),observed(indexKnown)',heterosig,heterol,heteroseps);   
                M1{ii,jj,kk}=reconstW(:);
                %[result2(ii,jj,kk),rs2(ii,jj,kk)]=GPR_core(observed2,dW,indexKnown2,length(indexKnown2),observed2(indexKnown2)',heterosig,heterol,heteroseps);            
                [result3(ii,jj,kk),rs3(ii,jj,kk),reconstW2]=GPR_core(observed3,dW,indexKnown3,length(indexKnown3),observed3(indexKnown3)',heterosig,heterol,heteroseps);     
                M2{ii,jj,kk}=reconstW2(:);
                %[result4(ii,jj,kk),rs4(ii,jj,kk)]=GPR_core(observed4,dW,indexKnown4,length(indexKnown4),observed4(indexKnown4)',heterosig,heterol,heteroseps);            
                %[result5(ii,jj,kk),rs5(ii,jj,kk)]=GPR_core(observed5,dW,indexKnown5,length(indexKnown5),observed5(indexKnown5)',heterosig,heterol,heteroseps);            
                
                
                istnoise=ndW_mean-dW;
                snr(ii,jj,kk)=norm(istnoise(:))/norm(dW(:));
                snr2(ii,jj,kk)=norm(istnoise(indexKnown))/norm(dW(indexKnown));
            end
    end
end

end
