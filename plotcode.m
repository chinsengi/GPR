close all
load('result');
result1 = squeeze(result(1, 2:4, 1:5, :));
result2 = squeeze(result(3, 2:4, 1:5, :));
result3 = squeeze(result(4, 2:4, 1:5, :));
load('result_nonaffine');
result1non = squeeze(result(1, 5, 1:5,:));
result2non = squeeze(result(3, 5, 1:5,:));
result3non = squeeze(result(4, 5, 1:5,:));
x = 10:10:50;
figure
gpPlot = [];
for i = 1:3
    data = squeeze(result1(i, :, :));
    y = mean(data,2);
    std_dev = std(data,0,2);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1', fliplr(curve2')];
    fill(x2, inBetween, 'b','facealpha',.05);
    hold on;
    gpPlot(i)= plot(x, y,'LineWidth', 2,'DisplayName', 'GP');
end
ylim([0.2*1e-3, 2.2*1e-3]); 
% plotgp(result1non, 'b');
legend(gpPlot, 'cd = 0.1', 'cd = 0.2', 'cd = 0.5');
title('paramset 1');
hold off
figure
for i = 1:3
    data = squeeze(result2(i, :, :));
    y = mean(data,2);
    std_dev = std(data,0,2);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1', fliplr(curve2')];
    fill(x2, inBetween, 'b','facealpha',.05);
    hold on;
    gpPlot(i)= plot(x, y,'LineWidth', 2,'DisplayName', 'GP');
end
ylim([0.5e-3, 7e-3]);
% plotgp(result2non,'b');
legend(gpPlot, 'cd = 0.1', 'cd = 0.2', 'cd = 0.5');
title('paramset 3');
hold off
figure
for i = 1:3
    data = squeeze(result3(i, :, :));
    y = mean(data,2);
    std_dev = std(data,0,2);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1', fliplr(curve2')];
    fill(x2, inBetween, 'b','facealpha',.05);
    hold on;
    gpPlot(i)= plot(x, y,'LineWidth', 2,'DisplayName', 'GP');
end
ylim([0.001, 0.02]);
% plotgp(result3non,'b');
legend(gpPlot,'cd = 0.1', 'cd = 0.2', 'cd = 0.5');
title('paramset 4');
hold off


function [] = plotgp(data, color)  
    x = 10:10:50;
    y = mean(data,2);
    std_dev = std(data,0,2);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1', fliplr(curve2')];
    fill(x2, inBetween, color,'facealpha',.05);
    hold on;
    gpPlot= plot(x, y,'LineWidth', 2,'DisplayName', 'GP');
    legend(gpPlot, 'nonaffine');
end