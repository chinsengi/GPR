% for i = filled+1:nent
%     cur = invindex(i);
%     [row, col] = unvec(cur, dim(1));
%     mu_pos(i-filled) = ret(row, col);
% end
% L = chol(nearestSPD(Kss))';
% alpha = L'\(L\mu_pos);
% mu_filled = Ks*alpha;
% for i = 1:filled
%     cur = invindex(i);
%     [row, col] = unvec(cur, dim(1));
%     ret(row, col) = mu_filled(i);
% end
% 
figure
paramno = 4;
load('./result/result');
load('./result/result_l');
load('./result/result_lsig');
result = squeeze(result(paramno,5,10:10:100, :));
result_l = squeeze(result_l(paramno, 5, 10:10:100, :));
result_lsig = squeeze(result_lsig(paramno, 5, 10:10:100, :));
% result(result>10) = 5e-4;
% result_l(result_l>10) = 5e-4;
% result_lsig(result_lsig>6) = 5e-4;
y = mean(result,2); % your mean vector;
x = 10:10:100;
std_dev = std(result,0,2);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1', fliplr(curve2')];
% fill(x2, inBetween, 'b','facealpha',.05);
hold on;
gpPlot= plot(x, y, 'b','LineWidth', 2,'DisplayName', 'GP');
y = mean(result_l,2); % your mean vector;
x = 10:10:100;
std_dev = std(result_l,0,2);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1', fliplr(curve2')];
% h = fill(x2, inBetween, 'r','facealpha',.05);
hold on;
gplPlot = plot(x, y, 'r', 'LineWidth', 2,'DisplayName', 'GP l');
y = mean(result_lsig,2); % your mean vector;
x = 10:10:100;
std_dev = std(result_lsig,0,2);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1', fliplr(curve2')];
% h = fill(x2, inBetween, 'g', 'facealpha',.05);
hold on;
% gplsigPlot = plot(x, y, 'g', 'LineWidth', 2,'DisplayName', 'GP l');
% legend([gpPlot, gplPlot, gplsigPlot], 'GP', 'GP l', 'GP lsig');
legend([gpPlot, gplPlot], 'GP', 'GP l');
hold off
