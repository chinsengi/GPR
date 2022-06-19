load('result_nmseE');
load('result_relerrE');
load('result_nmseI');
load('result_relerrI');
load('save_seed');
homepath='C:\Users\sengi\OneDrive - UW\lim\manuscript\figures';
font = 16;
save_flag = false;
for paramno = 1:1
    if paramno == 2
        continue
    end
    ntrials = 50;
    x = 10:10:100;
    savedir = 'C:\Users\Sengi\Desktop\lim\report\presentation';
    %plot relative error for excitatory synaptic plasticity 
    figure('Renderer', 'painters', 'Position', [10 10 900 900], 'DefaultAxesFontSize',font)
    data1 = squeeze(result_relerrE(1, paramno+1, 1:1, 1:10, 1:ntrials));
%     data1 = squeeze(result_nmseE(1, paramno+1, 1:1, 1:10, 1:ntrials));
    data2 = squeeze(result_relerrI(1, paramno+1, 1:1, 1:10, 1:ntrials));
%     data2 = squeeze(result_nmseI(1, paramno+1, 1:1, 1:10, 1:ntrials));
    hold on
    l1 = plot(x, mean(data1,2));
    l2 = plot(x, mean(data2, 2));
    plotmean(x,data1, l1.Color);
    plotmean(x, data2, l2.Color);
    yline(0.3730,'--', 'noise level Exc');
    yline(0.34, '--', 'noise level Inh');
    legend('Excitatory Synapse', 'Inhibitory Synapse');
    xlabel('# observations');
    ylabel('relative error');
%     title('Feed Forward EE Connection');
    if save_flag 
        saveas(gcf, [homepath, '/FF Double Perf.fig'])
    end
    
    % plot the relative performance at 50 samples
    figure('Renderer', 'painters', 'Position', [10 10 900 900], 'DefaultAxesFontSize',font)
    ratio = [2 .5 2/3 3/2 4 3 1 1/3 1/4];
    [sratio, ind] = sort(ratio);
    dataE = squeeze(result_relerrE([2:5 7:11], paramno+1, 1:1, 5, 1:ntrials));
    dataI = squeeze(result_relerrI([2:5 7:11], paramno+1, 1:1, 5, 1:ntrials));
    meanE = mean(dataE, 2);
    meanI = mean(dataI, 2);
    hold on
    lE = plot(sratio, meanE(ind));
    lI = plot(sratio, meanI(ind));
    plotmean(sratio, dataE(ind, :), lE.Color);
    plotmean(sratio, dataI(ind, :), lI.Color);
    set(gca, 'XScale', 'log')

    % plot the relative strength performance
    figure('Renderer', 'painters', 'Position', [10 10 900 900], 'DefaultAxesFontSize',font)
    data1 = squeeze(result_relerrE(2, paramno+1, 1:1, 1:10, 1:ntrials));
    data2 = squeeze(result_relerrE(3, paramno+1, 1:1, 1:10, 1:ntrials));
    data3 = squeeze(result_relerrE(4, paramno+1, 1:1, 1:10, 1:ntrials));
    data4 = squeeze(result_relerrE(5, paramno+1, 1:1, 1:10, 1:ntrials));
    data5 = squeeze(result_relerrE(7, paramno+1, 1:1, 1:10, 1:ntrials));
    data6 = squeeze(result_relerrE(8, paramno+1, 1:1, 1:10, 1:ntrials));
    data7 = squeeze(result_relerrE(9, paramno+1, 1:1, 1:10, 1:ntrials));
    data8 = squeeze(result_relerrE(10, paramno+1, 1:1, 1:10, 1:ntrials));
    data9 = squeeze(result_relerrE(11, paramno+1, 1:1, 1:10, 1:ntrials));
    data3(isnan(data3)) = 1;
    hold on
    l1 = plot(x, mean(data1, 2));
    l2 = plot(x, mean(data2, 2));
    l3 = plot(x, mean(data3, 2));
    l4 = plot(x, mean(data4, 2));
    l5 = plot(x, mean(data5, 2));
    l6 = plot(x, mean(data6, 2));
    l7 = plot(x, mean(data7, 2));
    l8 = plot(x, mean(data8, 2));
    l9 = plot(x, mean(data9, 2));
    plotmean(x, data1, l1.Color);
    plotmean(x, data2, l2.Color);
    plotmean(x, data3, l3.Color);
    plotmean(x, data4, l4.Color);
    plotmean(x, data5, l5.Color);
    plotmean(x, data6, l6.Color);
    plotmean(x, data7, l7.Color);
    plotmean(x, data8, l8.Color);
    plotmean(x, data9, l9.Color);
    yline(0.3730,'--', 'noise level Exc');
    yline(0.5089, '--', 'noise level Inh');
    legend('2:1', '1:2', '2:3', '3:2', '4:1', '3:1', '1:1', '1:3', '1:4');
    xlabel('# observations');
    ylabel('relative error');
%     title('Feed Forward EE Connection');
    if save_flag
        saveas(gcf, [homepath, '/ff double relative strength Excitatory.fig'])
    end
    
    figure('Renderer', 'painters', 'Position', [10 10 900 900], 'DefaultAxesFontSize',font)
    data1 = squeeze(result_relerrI(2, paramno+1, 1:1, 1:10, 1:ntrials));
    data2 = squeeze(result_relerrI(3, paramno+1, 1:1, 1:10, 1:ntrials));
    data3 = squeeze(result_relerrI(4, paramno+1, 1:1, 1:10, 1:ntrials));
    data4 = squeeze(result_relerrI(5, paramno+1, 1:1, 1:10, 1:ntrials));
    data5 = squeeze(result_relerrI(7, paramno+1, 1:1, 1:10, 1:ntrials));
    data6 = squeeze(result_relerrI(8, paramno+1, 1:1, 1:10, 1:ntrials));
    data7 = squeeze(result_relerrI(9, paramno+1, 1:1, 1:10, 1:ntrials));
    data8 = squeeze(result_relerrI(10, paramno+1, 1:1, 1:10, 1:ntrials));
    data9 = squeeze(result_relerrI(11, paramno+1, 1:1, 1:10, 1:ntrials));
    data3(isnan(data3)) = 1;
    hold on
    l1 = plot(x, mean(data1, 2));
    l2 = plot(x, mean(data2, 2));
    l3 = plot(x, mean(data3, 2));
    l4 = plot(x, mean(data4, 2));
    l5 = plot(x, mean(data5, 2));
    l6 = plot(x, mean(data6, 2));
    l7 = plot(x, mean(data7, 2));
    l8 = plot(x, mean(data8, 2));
    l9 = plot(x, mean(data9, 2));
    plotmean(x, data1, l1.Color);
    plotmean(x, data2, l2.Color);
    plotmean(x, data3, l3.Color);
    plotmean(x, data4, l4.Color);
    plotmean(x, data5, l5.Color);
    plotmean(x, data6, l6.Color);
    plotmean(x, data7, l7.Color);
    plotmean(x, data8, l8.Color);
    plotmean(x, data9, l9.Color);
    yline(0.3730,'--', 'noise level Exc');
    yline(0.5089, '--', 'noise level Inh');
    legend('2:1', '1:2', '2:3', '3:2', '4:1', '3:1', '1:1', '1:3', '1:4');
    xlabel('# observations');
    ylabel('relative error');
%     title('Feed Forward EE Connection');
    if save_flag
        saveas(gcf, [homepath, '/ff double relative strength Inhibitory.fig']);
    end

    % performance of EE network. 
    figure('Renderer', 'painters', 'Position', [10 10 900 900], 'DefaultAxesFontSize',font)
    data1 = squeeze(result_relerrE(6, paramno+1, 1:1, 1:10, 1:ntrials));
%     data1 = squeeze(result_nmseE(1, paramno+1, 1:1, 1:10, 1:ntrials));
    data2 = squeeze(result_relerrI(6, paramno+1, 1:1, 1:10, 1:ntrials));
%     data2 = squeeze(result_nmseI(1, paramno+1, 1:1, 1:10, 1:ntrials));
    hold on
    l1 = plot(x, mean(data1,2));
    l2 = plot(x, mean(data2, 2));
    plotmean(x,data1, l1.Color);
    plotmean(x, data2, l2.Color);
    yline(0.3730,'--', 'noise level Exc');
    yline(0.5089, '--', 'noise level Inh');
    legend('Excitatory Synapse', 'Excitatory Synapse 2');
    xlabel('# observations');
    ylabel('relative error');
%     title('Feed Forward EE Connection');
    if save_flag
        saveas(gcf, [homepath, '/ff double EE.fig']);
    end
    %plot R^2 for excitatory synaptic plasticity
%     figure('Renderer', 'painters', 'Position', [10 10 900 900])
%     data = squeeze(result_nmseE(1, paramno+1, 1:1, 1:10, 1:ntrials));
%     errorbar(x, [mean(data,2)], std(data,0,2));
%     hold on;
%     xlabel('# observations');
%     ylabel('NMSE');
%     title('Feed Forward EE Connection');
%     hold off
    
    %plot relative error for inhibitory synaptic plasticity 
    figure('Renderer', 'painters', 'Position', [10 10 900 900], 'DefaultAxesFontSize',font)
    data1 = squeeze(result_relerrI(1, paramno+1, 1:1, 1:10, 1:ntrials));
    hold on
    l1 = plot(x, mean(data1,2));
    plotmean(x,data1, l1.Color);
    legend('Inhibitory Synapse')
    xlabel('# observations');
    ylabel('relative error');
%     title('Feed Forward EI Connection');
    if save_flag
        saveas(gcf, [homepath, '/Feed Forward EI Connection.svg'])
    end
    
    %plot R^2 for inhibitory synaptic plasticity
%     figure('Renderer', 'painters', 'Position', [10 10 900 900])
%     data = squeeze(result_nmseI(1, paramno+1, 1:1, 1:10, 1:ntrials));
%     errorbar(x, [mean(data,2)], std(data, 0,2));
%     hold on;
%     xlabel('# observations');
%     ylabel('NMSE');
%     title('Feed Forward EI Connection');
%     hold off

    %plot median performance of excitatory map at 100 observation
    tmp = squeeze(result_relerrE(1, paramno+1, 1:1, 10, 1:ntrials-1));
    medianperf = median(tmp);
    wt = find(abs(tmp-medianperf)<1e-5, 1 ); %which trial has the median perfermance
    save_seed(1, paramno+1, wt)
    figure;
    imagesc('XData',irEy,'YData',irEx, 'CData',predE);
    cl=min(predE(:));
    ch=max(predE(:));
    cscale=normalize([cl ch/3 ch*2/3 ch],'range');
    cscale(1)=0;
    J=customcolormap(cscale, {'#f66e45','#ffffbb','#65c0ae','#5e4f9f'}, 64);
    colormap(J)
    colorbar
    caxis([cl ch])
    axis equal
    axis tight  
    xlabel('presynaptic firing rate (Hz)');
    ylabel('postsynaptic firing rate (Hz)')
    if save_flag
        saveas(gcf, [homepath, '/FFDpredE.fig'])
    end
    
    figure;
    imagesc('XData',irEy,'YData',irEx, 'CData',trueE);
    colormap(J)
    colorbar
    caxis([cl ch])
    axis equal
    axis tight  
    xlabel('presynaptic firing rate (Hz)');
    ylabel('postsynaptic firing rate (Hz)')
    if save_flag
        saveas(gcf, [homepath, '/FFDtrueE.fig'])
    end
%     saveas(gcf, [homepath, '/FFDtrueE.png'])

    % plot median performance of inhibitory map at 100 observations
    tmp = squeeze(result_relerrI(1, paramno+1, 1:1, 10, 1:ntrials-1));
    medianperf = median(tmp)
    wt = find(abs(tmp-medianperf)<1e-5, 1 ); %which trial has the median perfermance
    save_seed(1, paramno+1, wt)
    figure;
    imagesc('XData',irIy,'YData',irIx, 'CData',predI);
    cl=min(predE(:));
    ch=max(predE(:));
    cscale=normalize([cl ch/2 ch*3/4 ch],'range');
    cscale(1)=0; cscale(end) = 1;
    J=customcolormap(cscale, {'#f66e45','#ffffbb','#65c0ae','#5e4f9f'});
    colormap(J)
    colorbar
    caxis([cl ch])
    axis equal
    axis tight  
    xlabel('presynaptic firing rate (Hz)');
    ylabel('postsynaptic firing rate (Hz)')
    if save_flag
        saveas(gcf, [homepath, '/FFDpredI.fig'])
    end
%     saveas(gcf, [homepath, '/FFDpredI.png'])
    
    figure;
    imagesc('XData',irIy,'YData',irIx, 'CData',trueI);
    colormap(J)
    colorbar
    caxis([cl ch])
    axis equal
    axis tight  
    xlabel('presynaptic firing rate (Hz)');
    ylabel('postsynaptic firing rate (Hz)')
    if save_flag
        saveas(gcf, [homepath, '/FFDtrueI.fig'])
    end
%     saveas(gcf, [homepath, '/FFDI.png'])
end

hold off
x = 0:5:(max([Rate_Novel; Rate_Fam])+5);
figure(1);
y = histogram(Rate_Novel, x);
figure(2);
stairs(x(1:end-1),y.Values);
hold on
figure(1);
y = histogram(Rate_Fam, x);
figure(2);
stairs(x(1:end-1),y.Values);
legend('Rate Novel', 'Rate Familiar');
ylim([0,25])
% legend(gpPlot, '1000 neurons', '2000 neurons', '5000 neurons');