function [dW,indexKnown]=voronoisampling(dim,obsinit,ndW)
dW = zeros(dim); % dW with missing value
try
    coordvoro=floor(mycvt_2d_sampling(obsinit,20,2500)*50);
catch
    coordvoro=floor(mycvt_2d_sampling(obsinit,20,2500)*50);
end
coordvoro(coordvoro<=0)=1;
dWt=zeros(size(ndW));
indt=sub2ind(dim,coordvoro(:,2),coordvoro(:,1));
dWt(indt)=ndW(indt);
rowsselect=coordvoro(:,1);
init=coordvoro(:,2);
[rowsselect,init]=find(dWt~=0);    
for j = 1:obsinit
    i=rowsselect(j);
    dW(i,init(j)) = ndW(i,init(j));
end
indexKnown=find(dW~=0);