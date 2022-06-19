load('plasticity map.mat')
load('Novel_Fam.mat')

Rate_Novel = tmp(:,1);
Rate_Fam = tmp(:,2);
InputChange = Rate_Fam-Rate_Novel;

[U,S,V] = svd(X);
% U = (1:70)';
vq = interp1(1:70,U(:,1),Rate_Novel,'spline');

figure(1)
plot(U(:,1));hold on;
plot(Rate_Novel,vq,'.');hold off

coeff = sum(InputChange.*vq)/sum(vq.*vq);

figure(2)
plot(Rate_Novel,InputChange,'o');hold on
plot(1:70,coeff*U(:,1));hold off