load('result_reg_rel3');
load('result_reg_rel4');
load('result_relerr')
load('save_seed');
nnlist = [200, 500, 1000];
font = 16;
homepath='C:\Users\sengi\OneDrive - UW\lim\manuscript\figures';
for paramno = 1:1
    if paramno == 2
        continue
    end
    ntrials = 24;
    %plot relative error
    figure('Renderer', 'painters', 'Position', [10 10 700 700], 'DefaultAxesFontSize',font)
    ob_range = 1:10;
    x = ob_range*10;
    data1 = squeeze(result_relerr(1, paramno+1, 1, ob_range, 1:ntrials));
    data2 = squeeze(result_reg_rel3(2, paramno+1, ob_range, 1:ntrials));
    data3 = squeeze(result_relerr(3, paramno+1, 1, ob_range, 1:ntrials));
    data4 = squeeze(result_relerr(5, paramno+1, 1, ob_range, 1:ntrials));
    data5 = squeeze(result_relerr(2, paramno+1, 1, ob_range, 1:ntrials));
    data6 = squeeze(result_relerr(4, paramno+1, 1, ob_range, 1:ntrials));
    data7 = squeeze(result_reg_rel3(1, paramno+1, ob_range, 1:ntrials));
    hold on
    l1 = plot(x, mean(data1,2));
    l2 = plot(x, mean(data2,2));
    l3 = plot(x, mean(data3,2));
    l4 = plot(x, mean(data4,2));
    l5 = plot(x, mean(data5,2));
%     l6 = plot(x, mean(data6,2)); 
    l7 = plot(x, mean(data7,2));
    plotmean(x,data1, l1.Color);
    plotmean(x,data2, l2.Color);
    plotmean(x,data3, l3.Color);
    plotmean(x,data4, l4.Color);
    plotmean(x,data5, l5.Color);
    plotmean(x, data7, l7.Color)
%     plotmean(x,data6, l6.Color);
%     breakyaxis([1 10]);
    ylim([0,0.5]);
    yline(0.3863, '--', 'noise level');
    legend('GPR', 'TE3 missing 10% connections', 'input noise \sigma = 0.1',...
        'input noise \sigma = 0.5', 'missing 10% connections', 'TE3', 'noise level');
%     legend('GPR', 'active GPR')
    xlabel('number of observations');
    ylabel('relative error');
%     title('FF-E');
    hold off
    xlim([10, 100]);
    saveas(gcf, [homepath, '/FF-E.png'])
    saveas(gcf, [homepath, '/FF-E.fig'])

    % plot R^2 
%     data = squeeze(result_rsq(1, paramno+1, 1:1, 1:10, 1:ntrials));
%     figure('Renderer', 'painters', 'Position', [10 10 900 900]);
%     x = 10:10:100;
%     errorbar(x, mean(data, 2), std(data,0,2));
%     hold on
%     xlabel('# observations');
%     ylabel('R^2');
%     title('FF-E');
%     hold off
    
    % plot input noise
%     data1 = squeeze(result_relerr(3, paramno+1, 1, 1:10, 1:ntrials));
%     data2 = squeeze(result_relerr(4, paramno+1, 1, 1:10, 1:ntrials));
%     data3 = squeeze(result_relerr(1, paramno+1, 1, 1:10, 1:ntrials));
%     figure('Renderer', 'painters', 'Position', [10 10 900 900]);
%     x = 10:10:100;
%     hold on; 
%     boxchart(data1');
%     boxchart(data2');
%     boxchart(data3');
%     xticklabels(string(x));
%     xlabel('# observations');
%     ylabel('relative error');
% %     title('FF-E input noise \sigma = 0.1');
%     legend('input noise w/o correction', 'input noise with correction', 'no input noise');
%     hold off
%     saveas(gcf, [homepath, '/FF-E input noise sigma 0.1.svg'])
%     
%     data1 = squeeze(result_relerr(5, paramno+1, 1, 1:10, 1:ntrials));
%     data2 = squeeze(result_relerr(6, paramno+1, 1, 1:10, 1:ntrials));
%     figure('Renderer', 'painters', 'Position', [10 10 900 900]);
%     x = 10:10:100;
%     hold on; 
%     boxchart(data1');
%     boxchart(data2');
%     boxchart(data3');
%     xticklabels(string(x));
%     xlabel('# observations');
%     ylabel('relative error');  
% %     title('FF-E input noise \sigma = 0.5');
%     legend('input noise w/o correction', 'input noise with correction', 'no input noise');
%     hold off
%     saveas(gcf, [homepath, '/FF-E input noise sigma 0.5.svg'])
    
    % effect of number of projections on performance 
    data1 = squeeze(result_relerr(9, paramno+1, 1, 1:10, 1:ntrials));
    data2 = squeeze(result_relerr(10, paramno+1, 1, 1:10, 1:ntrials));
    data3 = squeeze(result_relerr(11, paramno+1, 1, 1:10, 1:ntrials));
    data4 = squeeze(result_relerr(12, paramno+1, 1, 1:10, 1:ntrials));
    data5 = squeeze(result_relerr(1, paramno+1, 1, 1:10, 1:ntrials));
    figure('Renderer', 'painters', 'Position', [10 10 900 900], 'DefaultAxesFontSize',font)
    x = 10:10:100;
    hold on; 
    l1 = plot(x, mean(data1,2));
    l2 = plot(x, mean(data2,2));
    l3 = plot(x, mean(data3,2));
    l4 = plot(x, mean(data4,2));
    l5 = plot(x, mean(data5,2));
    plotmean(x,data1, l1.Color);
    plotmean(x,data2, l2.Color);
    plotmean(x,data3, l3.Color);
    plotmean(x,data4, l4.Color);
    plotmean(x,data5, l5.Color);
    xlabel('# observations');
    ylabel('relative error');
    yline(0.3863, '--', 'noise level');
%     title('FF-E input noise \sigma = 0.1');
    legend('100 pre neurons', '200 pre neurons', '300 pre neurons', '400 pre neurons', '500 pre neurons', 'noise level');
    gca.FontSize= font;
    hold off
    xlim([10, 100]);
    saveas(gcf, [homepath, '/FF-E projection.png'])
    saveas(gcf, [homepath, '/FF-E projection.fig'])
    
    % plot missing connection 
%     data1 = squeeze(result_relerr(1, paramno+1, 1, 1:10, 1:ntrials));
%     data2 = squeeze(result_relerr(2, paramno+1, 1, 1:10, 1:ntrials));
%     data3 = squeeze(result_relerr(7, paramno+1, 1, 1:10, 1:ntrials));
%     data4 = squeeze(result_reg_rel4(paramno+1, 1:10, 1:ntrials));
%     figure('Renderer', 'painters', 'Position', [10 10 900 900], 'DefaultAxesFontSize',font)
%     x = 10:10:100;
%     hold on
%     l1 = plot(x, mean(data1,2));
%     l2 = plot(x, mean(data2,2));
%     l3 = plot(x, mean(data3,2));
%     l4 = plot(x, mean(data4,2));
%     plotmean(x,data1, l1.Color);
%     plotmean(x,data2,l2.Color);
%     plotmean(x,data3, l3.Color);
%     plotmean(x, data4, l4.Color);
%     xlabel('# observations');
%     ylabel('relative error');
%     legend('Full', '10% missing connection', '20% missing connection', 'TE3');
% %     title('FF-E missing connection');
%     gca.FontSize= font;
%     hold off
%     saveas(gcf, [homepath, '/FF-E missing connection.svg'])
%     saveas(gcf, [homepath, '/FF-E missing connection.png'])

    %plot median performance at 50 observation
    expno = 1;
    model_selection_method = 1;
    tmp = squeeze(result_relerr(expno, paramno+1, model_selection_method,...
                        5, 1:ntrials-1));
    medianperf = median(tmp);
    wt = find(abs(tmp-medianperf)<1e-5, 1); %which trial has the median perfermance
    save_seed(expno, paramno+1, wt)
    figure;
    imagesc('XData',iry,'YData',irx, 'CData',ret);
    cl=min(ret(:));
    ch=max(ret(:));
    cscale=normalize([cl ch/3 ch*2/3 ch],'range');
    cscale(1)=0;
    J=customcolormap(cscale, {'#f66e45','#ffffbb','#65c0ae','#5e4f9f'});
    colormap(J)
    colorbar
    caxis([cl ch])
    axis equal
    axis tight  
    xlabel('r_{pre} (Hz)');
    ylabel('r_{post} (Hz)')
    saveas(gcf, [homepath, '/FFSpredE.fig'])
%     saveas(gcf, [homepath, '/FFSpredE.png'])
%     
    figure;
    imagesc('XData',iry,'YData',irx, 'CData',Xp);
    colormap(J)
    colorbar
    caxis([cl ch])
    axis equal
    axis tight  
    xlabel('r_{pre} (Hz)');
    ylabel('r_{post} (Hz)')
    saveas(gcf, [homepath, '/FFStrueE.fig'])
%     saveas(gcf, [homepath, '/FFStrueE.png'])
    
end

figure;
fam_hist = histogram(Rate_Fam, linspace(0, 80, 20));
fam_edge = fam_hist.BinEdges;
fam_value = fam_hist.Values;
nov_hist = histogram(Rate_Novel, linspace(0, 80, 20));
nov_edge = nov_hist.BinEdges;
nov_value = nov_hist.Values;
figure('Renderer', 'painters', 'Position', [10 10 900 900], 'DefaultAxesFontSize',40)
stairs(fam_edge, [fam_value, 0]);
hold on; 
stairs(nov_edge, [nov_value, 0]);
legend('Rate Familiar', 'Rate Novel');
saveas(gcf, [homepath, '/FF-E_rate_dist.fig'])