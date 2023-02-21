% load('./result/result_rel');
% load('./result/result_rel_reg3')
% load('./result/result_rel_reg4')
for paramno = 1:1
    if paramno == 2
        continue
    end
    x = 10:10:100;
    result_1 = squeeze(result(1, paramno, 2, x, 1:50));
    result_2 = squeeze(result(2, paramno, 2, x, 1:50));
%     result_1 = squeeze(result_logp(1,paramno,1, x, 1:50));
%     result_2 = squeeze(result_logp(2,paramno,1, x, 1:50));

    savedir = 'C:\Users\Sengi\Desktop\lim\report\presentation';
    %plot relative error
    method = 'gpr';
    figure('Renderer', 'painters', 'Position', [10 10 900 900])
    data = mean(result_1,2);
    data_fixed = mean(result_2,2);
    plot(x, [data, data_fixed]);
    
    hold on;
    title(sprintf(strcat(method,' nonaffine param %d'), paramno));
    xlabel('# observations');
    ylabel('relative error');
    fname = strcat(method, '_nonaffine_p', string(paramno), '_rel');
%         ylevels = yticks();
%         if ylevels(2) - ylevels(1)> 0.1
%          yticks(ylevels(1):0.1:ylevels(end));
%         end
    saveas(gcf, fullfile(savedir, fname), 'png');
    hold off

%     method = 'reg3';
%     figure('Renderer', 'painters', 'Position', [10 10 900 1200])
%     data = result_2;
%     boxplot(data', x);
%     hold on;
%     title(sprintf(strcat(method,' nonaffine param %d'), paramno));
%     xlabel('# observations');
%     ylabel('relative error');
%     fname = strcat(method,'_nonaffine_p', string(paramno), '_rel');
%         ylevels = yticks();
%         if ylevels(2) - ylevels(1)> 0.1
%          yticks(ylevels(1):0.1:ylevels(end));
%         end
%     saveas(gcf, fullfile(savedir, fname), 'png');
%     hold off
%     
%     method = 'reg4';
%     figure('Renderer', 'painters', 'Position', [10 10 900 1200])
%     data = result_3;
%     boxplot(data', x);
%     hold on;
%     title(sprintf(strcat(method,' nonaffine param %d'), paramno));
%     xlabel('# observations');
%     ylabel('relative error');
%     fname = strcat(method,'_nonaffine_p', string(paramno), '_rel');
%         ylevels = yticks();
%         if ylevels(2) - ylevels(1)> 0.1
%          yticks(ylevels(1):0.1:ylevels(end));
%         end
    hold off
end
% legend(gpPlot, '1000 neurons', '2000 neurons', '5000 neurons');