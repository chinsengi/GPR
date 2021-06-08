c = 0;
% for i = 1:500
a = (2*rand(10000,1)+4)/500;
b = exprnd(0.619, 100000, 1);
c = a.*Diff_Rate;
% end
histfit(c, 100, 'gamma')