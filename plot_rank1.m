% load('result_reg_rel3');
% load('result_reg_rel4');
% load('result_rel')
nnlist = [200, 500, 1000];
for paramno = 1:1
    if paramno == 2
        continue
    end
    nsamples = 50;
    result_1 = squeeze(result_rel_cval(paramno+1, 1:2, 2:10, 1:nsamples));
%     result_1 = squeeze(result_rel_input_vanilla(paramno+1, 1:2, 2:10, 1:nsamples));
%     result_1 = squeeze(result_rel(paramno+1, 1:2, 2:10, 1:nsamples));
%     result_2 = squeeze(result_reg_rel3(paramno+1,1:2, 2:10, 1:nsamples));
%     result_3 = squeeze(result_reg_rel4(paramno+1,1:2, 2:10, 1:nsamples));
%     result_2 = squeeze(mean_diff(paramno+1, 1:2, 1:15, :));
    x = 20:10:100;
    savedir = 'C:\Users\Sengi\Desktop\lim\report\presentation';
    for i = 2:2
        %plot relative error
        method = 'gpr extra term with MAP';
        figure('Renderer', 'painters', 'Position', [10 10 900 1200])
        data = squeeze(result_1(i, :, :));
        boxplot(data', x);
        hold on;
        title(sprintf(strcat(method,'%d neurons param %d'), nnlist(i), paramno));
        xlabel('# observations');
        ylabel('relative error');
        fname = strcat(method, '_', string(nnlist(i)), 'nn_p', string(paramno), '_rel');
%         ylevels = yticks();
%         if ylevels(2) - ylevels(1)> 0.1
%          yticks(ylevels(1):0.1:ylevels(end));
%         end
        saveas(gcf, fullfile(savedir, fname), 'png');
        hold off
        
%         method = 'reg3';
%         figure('Renderer', 'painters', 'Position', [10 10 900 1200])
%         data = squeeze(result_2(i, :, :));
%         boxplot(data', x);
%         hold on;
%         title(sprintf(strcat(method,' %d neurons param %d'), nnlist(i), paramno));
%         xlabel('# observations');
%         ylabel('relative error');
%         fname = strcat(method, '_', string(nnlist(i)), 'nn_p', string(paramno), '_rel');
%         ylevels = yticks();
%         if ylevels(2) - ylevels(1)> 0.1
%          yticks(ylevels(1):0.1:ylevels(end));
%         end
        saveas(gcf, fullfile(savedir, fname), 'png');
        hold off
        
%         method = 'reg4';
%         figure('Renderer', 'painters', 'Position', [10 10 900 1200])
%         data = squeeze(result_3(i, :, :));
%         boxplot(data', x);
%         hold on;
%         title(sprintf(strcat(method,' %d neurons param %d'), nnlist(i), paramno));
%         xlabel('# observations');
%         ylabel('relative error');
%         fname = strcat(method, '_', string(nnlist(i)), 'nn_p', string(paramno), '_rel');
%         ylevels = yticks();
%         if ylevels(2) - ylevels(1)> 0.1
%          yticks(ylevels(1):0.1:ylevels(end));
%         end
        saveas(gcf, fullfile(savedir, fname), 'png');
        hold off
        %plot mean difference
%         figure
%         data = squeeze(result_2(i, :, :));
%         boxplot(data', x);
%         hold on;
%         title(sprintf('%d neurons param %d', nnlist(i), paramno));
%         xlabel('# observations');
%         ylabel('mean difference');
%         fname = strcat(method, '_', string(nnlist(i)), 'nn_p', string(paramno), '_md');
%         saveas(gcf, fullfile(savedir, fname), 'png');
    %     ylabel('relative error');
    %     y = mean(data,2);
    %     std_dev = std(data,0,2);
    %     curve1 = y + std_dev;
    %     curve2 = y - std_dev;
    %     x2 = [x, fliplr(x)];
    %     inBetween = [curve1', fliplr(curve2')];
    % %     fill(x2, inBetween, 'b','facealpha',.05);
    %     hold on;
    %     gpPlot(i)= plot(x, y,'LineWidth', 2,'DisplayName', 'GP');
    end
end
% legend(gpPlot, '1000 neurons', '2000 neurons', '5000 neurons');