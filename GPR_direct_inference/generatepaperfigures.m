clear;
addpath('D:\Documents\Lim_Lab\apg');
addpath('D:\Documents\Lim_Lab\baby_project_code');
addpath('D:\Documents\Lim_Lab\baby_project_code\utils');
addpath('D:\Documents\Lim_Lab\LMaFit-adp');
addpath('E:\Documents\Lim_Lab\MatrixCompletion');
addpath('E:\Documents\Lim_Lab\customcolormap');
data2=load('D:\Documents\Lim_Lab\CalciumModel\dWmatrices2\CpreCpost1.mat');
addpath('E:\Documents\Lim_Lab\DRM_GPR');
data2=load('CpreCpost1.mat');
dWs=data2.dWs;
dWvars=data2.dWvars;
f2=1:50;
data3=load('paramset1.mat');
dWs{8}=data3.dW1(1:2:end,1:2:end);
dWvars{8}=data3.dWvar(1:2:end,1:2:end);
homepath='E:\Documents\Lim_Lab\hetero-gaussian-regression2\finalfigures';
%% Figure 1
result=cell(length(dWs),1);
result2=cell(length(dWs),1);
result3=cell(length(dWs),1);
result4=cell(length(dWs),1);
rs=cell(length(dWs),1);
rs2=cell(length(dWs),1);
rs3=cell(length(dWs),1);
rs4=cell(length(dWs),1);
snr=cell(length(dWs),1);
snr2=cell(length(dWs),1);
nsim=1;
ntrial=1;
srs=20;
M1s=cell(length(dWs),1);
M2s=cell(length(dWs),1);
for mm=8
    X=dWs{mm,1};
    dWvar=dWvars{mm,1};
    [result{mm},result2{mm},result3{mm},rs{mm},rs2{mm},rs3{mm},snr{mm},snr2{mm},M1s{mm},M2s{mm}]=GPR_TE_DRM(X,dWvar,false,false,false,nsim,ntrial,srs);
    %[result{mm},result2{mm},result3{mm},result4{mm},rs{mm},rs2{mm},rs3{mm},rs4{mm},snr{mm},snr2{mm},M1s{mm},M2s{mm}]=GPR_TE_DRM_SVT(X,dWvar,false,false,false,nsim,ntrial,srs);
end
dW=(dWs{8}-0.5)/0.5;
figure;
tiledlayout(1,4)
nexttile;
imagesc('XData',f2,'YData',f2, 'CData',dW);
J = customcolormap_preset('pasteljet');
%J=customcolormap(linspace(0,1,5),{'#FF0000','#FFFF00','#00FF00','#00FFFF','#0000FF'});
%J=customcolormap([0 0.25 0.5 0.75 1],{'#FF0000','#FFFF00','#00FF00','#00FFFF','#0000FF'});

cl=min(dW(:));
ch=max(dW(:));
%% Figure1A (breakpoint in GRP_TE)
%dW=(dWs{8}-0.5)/0.5;
figure;
tiledlayout(1,4)
nexttile;
imagesc('XData',f2,'YData',f2, 'CData',dW);
J = customcolormap_preset('pasteljet');
%J=customcolormap(linspace(0,1,5),{'#FF0000','#FFFF00','#00FF00','#00FFFF','#0000FF'});
%J=customcolormap([0 0.25 0.5 0.75 1],{'#FF0000','#FFFF00','#00FF00','#00FFFF','#0000FF'});
cl=min(dW(:));
ch=max(dW(:));
cscale=normalize([cl  0 ch/3 ch*2/3 ch],'range');
cscale(1)=0;
J=customcolormap(cscale, {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'});
colormap(J)
box off;
axis equal
axis tight
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')

nexttile
imagesc('XData',f2,'YData',f2, 'CData',tdW);
caxis([cl ch]) 
%colormap(J)
axis equal
axis tight
box off;
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')

observed(observed==0)=NaN;
imAlpha=ones(size(observed));
imAlpha(isnan(observed))=0;
nexttile
imagesc('XData',f2,'YData',f2, 'CData',observed,'AlphaData',imAlpha);
caxis([cl ch]) 
colormap(J)
box off;
axis equal
axis tight
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')

nexttile
imagesc('XData',f2,'YData',f2, 'CData',reconstW);
caxis([cl ch]) 
colormap(J)
box off;
axis equal
axis tight
colorbar
xlabel('presynaptic firing rate (Hz)')
ylabel('postsynaptic firing rate (Hz)')

set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
%% Figure 1B
result=cell(length(dWs),1);
result2=cell(length(dWs),1);
result3=cell(length(dWs),1);
result4=cell(length(dWs),1);
rs=cell(length(dWs),1);
rs2=cell(length(dWs),1);
rs3=cell(length(dWs),1);
rs4=cell(length(dWs),1);
snr=cell(length(dWs),1);
snr2=cell(length(dWs),1);
nsim=50;
ntrial=1;
srs=10:10:100;
M1s=cell(length(dWs),1);
M2s=cell(length(dWs),1);
M_samples=cell(length(dWs),1);
for mm=8
    X=dWs{mm,1};
    dWvar=dWvars{mm,1};
    [result{mm},result2{mm},result3{mm},rs{mm},rs2{mm},rs3{mm},snr{mm},snr2{mm},M1s{mm},M2s{mm}]=GPR_TE_DRM(X,dWvar,false,false,false,nsim,ntrial,srs);
    %[result{mm},result2{mm},result3{mm},result4{mm},rs{mm},rs2{mm},rs3{mm},rs4{mm},snr{mm},snr2{mm},M1s{mm},M2s{mm},M_samples{mm}]=GPR_TE_DRM_SVT(X,dWvar,false,false,false,nsim,ntrial,srs);
end
save(fullfile(homepath,'fig1b.mat'),'srs','ntrial','nsim','result','result2','rs','rs2','snr','snr2','M1s','M2s');
%%
figure;

dy=0.4;
for mm=8
    mean_result=mean(squeeze(result{mm}),2);
    mean_result2=mean(squeeze(result2{mm}),2);
    mean_result3=mean(squeeze(result3{mm}),2);
    mean_result4=mean(squeeze(result4{mm}),2);
    mean_snr=mean(squeeze(snr{mm}),2);
    mean_snr2=mean(squeeze(snr2{mm}),2);
    hold on;
    l1=plot(srs,mean_result);
    l2=plot(srs,mean_result2);
    l3=plot(srs,mean_result3);
    l4=plot(srs,mean_result4-dy);
    plotmean(srs,squeeze(result{mm}),l1.Color);
    plotmean(srs,squeeze(result2{mm}),l2.Color);
    plotmean(srs,squeeze(result3{mm}),l3.Color);
    plotmean(srs,squeeze(result4{mm}-dy),l4.Color);
    plot(srs,mean_snr,'--k');
    %yline(0.15,'--')
    %plotmean(srs,squeeze(result3{mm}),l3.Color);

end
plot(50,mean_result(5),'*','Color',l1.Color)
plot(50,mean_result2(5),'*','Color',l2.Color)
xlabel('number of observed points');
ylabel('relative error');
axis tight
lgd=legend('GPR','TE','DRM','SVT');
lgd.FontSize = 8;
box off
x0=100;
y0=200;
width=130.2;
height=102.8;
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
ylim([0.1 0.6])
yticks([0.15:0.05:0.4, 0.5:0.05:0.6])
yticklabels([0.15:0.05:0.4,(0.5:0.05:0.6)+dy])
axes('Position',[.1 .65 .05 .05]);
px=[1 5];
py1=[1 2];
height=1;
py2=py1+height;
plot(px,py1,'k','LineWidth',2);hold all;
plot(px,py2,'k','LineWidth',2);hold all;
fill([px flip(px)],[py1 flip(py2)],'w','EdgeColor','none');
box off;
axis off;
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
saveas(gcf,fullfile(homepath,'fig1b.eps'))
saveas(gcf,fullfile(homepath,'fig1b.png'))
savefig(fullfile(homepath,'fig1b.fig'))
%%
figure;
for mm=8
    mean_result=mean(squeeze(result{mm}),2);
    mean_result2=mean(squeeze(result2{mm}),2);
    mean_result3=mean(squeeze(result3{mm}),2);
    mean_snr=mean(squeeze(snr{mm}),2);
    mean_snr2=mean(squeeze(snr2{mm}),2);
    nexttile
    hold on;
    boxchart(squeeze(result{mm})','BoxFaceAlpha',0.2,'BoxWidth',0.4,'LineWidth',1,'MarkerSize',5)
    boxchart(squeeze(result2{mm})','BoxFaceAlpha',0.2,'BoxWidth',0.4,'LineWidth',1,'MarkerSize',5)
    plot(srs/10,mean_snr,'--k');
    axis tight
end
box off
xticklabels(10:10:100);
xlabel('number of observed points');
ylabel('relative error');
legend('GPR','TE')
x0=100;
y0=200;
width=130.2;
height=102.8;

%set(gcf,'units','points','position',[x0,y0,width,height])
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
saveas(gcf,fullfile(homepath,'fig1b.eps'))
saveas(gcf,fullfile(homepath,'fig1b.png'))
savefig(fullfile(homepath,'fig1b.fig'))
%%
figure;
tiledlayout(2,1)
for mm=8
    rtemp=squeeze(result{mm});
    [~,medianindex]=min(abs(rtemp(2,:)-median(rtemp(2,:))));
end
imgs=M1s{8};
nexttile;
imagesc('XData',f2,'YData',f2, 'CData',reshape(imgs{1,2,medianindex},[50,50]));
caxis([cl ch]) 
colormap(J)
axis equal
axis tight

for mm=8
    rtemp=squeeze(result2{mm});
    [~,medianindex]=min(abs(rtemp(2,:)-median(rtemp(2,:))))
end
imgs=M2s{8};
nexttile
imagesc('XData',f2,'YData',f2, 'CData',reshape(imgs{1,2,medianindex},[50,50]));
caxis([cl ch]) 
colormap(J)
axis equal
axis tight
saveas(gcf,fullfile(homepath,'fig1c.eps'))
saveas(gcf,fullfile(homepath,'fig1c.png'))
savefig(fullfile(homepath,'fig1c.fig'))
%%
figure;
for mm=8
    rtemp=squeeze(result2{mm});
    [~,medianindex]=min(abs(rtemp(2,:)-median(rtemp(2,:))));
    imgs=M2s{mm};
    a3=nexttile;
    MMM=M_samples{mm};
    observed=MMM{1,2,medianindex};
    
    observed(observed==0)=NaN;
    imAlpha=ones(size(observed));
    imAlpha(isnan(observed))=0;
    imagesc('XData',f2,'YData',f2, 'CData',observed,'AlphaData',imAlpha);

    %title(num2str(rtemp(10,medianindex)));
    caxis([cl ch]) 
    axis equal
    axis tight
    colormap(a3,J)
    error1(mm)=rtemp(10,medianindex);
    
end

set(gcf,'units','points','position',[x0,y0,width,height])
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
saveas(gcf,fullfile(homepath,'newfig2_voronoisample.eps'))
saveas(gcf,fullfile(homepath,'newfig2_voronoisample.png'))
savefig(fullfile(homepath,'newfig2_voronoisample.fig'))
%%
figure;
hold on;
for mm=8
    mean_result=mean(squeeze(result{mm}),2);
    mean_result2=mean(squeeze(result2{mm}),2);
    mean_snr=mean(squeeze(snr{mm}),2);
    mean_snr2=mean(squeeze(snr2{mm}),2);
    %scatter(result{mm}(:),result2{mm}(:),8,cc(:,mm)','filled');
    scatter(result{mm}(:),result2{mm}(:),8,snr2{mm}(:),'filled');
end
plot([0,1.2],[0,1.2],'--k')
xlabel('GPR relative error')
ylabel('TE relative error')
axis tight
% caxis([0.1 0.8])
h=colorbar;%('XTick', 0.1:0.2:0.8);
colormap(gray*0.8)
ylim([0,1])
ylabel(h, 'relative noise level')
x0=200;
y0=200;
width=130.2;
height=102.8;

set(gcf,'units','points','position',[x0,y0,width,height])
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
box off
saveas(gcf,fullfile(homepath,'newfig2_voronoi.eps'))
saveas(gcf,fullfile(homepath,'newfig2_voronoi.png'))
savefig(fullfile(homepath,'newfig2_voronoi.fig'))
%%
figure;
for mm=8
    rtemp=squeeze(result2{mm});
    [~,medianindex]=min(abs(rtemp(2,:)-median(rtemp(2,:))));
    imgs=M2s{mm};
    a3=nexttile;

    M_r=reshape(squeeze(imgs{1,2,medianindex}),[50 50]);
    imagesc('XData',f2,'YData',f2, 'CData',M_r);

    %title(num2str(rtemp(10,medianindex)));
    caxis([cl ch]) 
    axis equal
    axis tight
    colormap(a3,J)
    error1(mm)=rtemp(10,medianindex);
    
end

set(gcf,'units','points','position',[x0,y0,width,height])
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
saveas(gcf,fullfile(homepath,'newfig2_voronoirecovery.eps'))
saveas(gcf,fullfile(homepath,'newfig2_voronoirecovery.png'))
savefig(fullfile(homepath,'newfig2_voronoirecovery.fig'))

%% Figure 1d
result=cell(length(dWs),1);
result2=cell(length(dWs),1);
result3=cell(length(dWs),1);
rs=cell(length(dWs),1);
rs2=cell(length(dWs),1);
rs3=cell(length(dWs),1);
snr=cell(length(dWs),1);
snr2=cell(length(dWs),1);
nsim=50;
ntrial=20;
srs=20;
M1s=cell(length(dWs),1);
M2s=cell(length(dWs),1);
for mm=8
    X=dWs{mm,1};
    dWvar=dWvars{mm,1};
    [result{mm},result2{mm},result3{mm},rs{mm},rs2{mm},rs3{mm},snr{mm},snr2{mm},M1s{mm},M2s{mm}]=GPR_TE(X,dWvar,false,false,false,nsim,ntrial,srs);
end
save(fullfile(homepath,'fig1d'),'srs','ntrial','nsim','result','result2','rs','rs2','snr','snr2','M1s','M2s');
%%
trials=1:ntrial;

figure;

for mm=8
    mean_result=mean(squeeze(result{mm}),2);
    mean_result2=mean(squeeze(result2{mm}),2);
    mean_result3=mean(squeeze(result3{mm}),2);
    mean_snr=mean(squeeze(snr{mm}),2);
    mean_snr2=mean(squeeze(snr2{mm}),2);
    nexttile([1,2])
    hold on;
%     l1=plot(trials,mean_result);
%     l2=plot(trials,mean_result2);
%    % l3=plot(trials,mean_result3);
%     plotmean(trials,squeeze(result{mm}),l1.Color);
%     plotmean(trials,squeeze(result2{mm}),l2.Color);
    boxchart(squeeze(result{mm})')
    boxchart(squeeze(result2{mm})')
    plot(trials,mean_snr,'--k');
    %plotmean(trials,squeeze(result3{mm}),l3.Color);
    axis tight
end

xlabel('number of trials');
ylabel('relative error');
legend('GPR','TE','noise level')
x0=100;
y0=200;
width=130.2;
height=102.8;
box off;

set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')

% for mm=8
%     rtemp=squeeze(result{mm});
%     [~,medianindex]=min(abs(rtemp(10,:)-median(rtemp(10,:))));
% end
% imgs=M1s{8};
% nexttile;
% imagesc('XData',f2,'YData',f2, 'CData',reshape(imgs{10,1,medianindex},[50,50]));
% caxis([cl ch]) 
% colormap(J)
% axis equal
% axis tight
% 
% for mm=8
%     rtemp=squeeze(result2{mm});
%     [~,medianindex]=min(abs(rtemp(10,:)-median(rtemp(10,:))));
% end
% imgs=M2s{8};
% nexttile;
% imagesc('XData',f2,'YData',f2, 'CData',reshape(imgs{10,1,medianindex},[50,50]));
% caxis([cl ch]) 
% colormap(J)
% box off;
% axis equal
% axis tight
% colorbar
% xlabel('presynaptic firing rate (Hz)')
% ylabel('postsynaptic firing rate (Hz)')

saveas(gcf,fullfile(homepath,'newfigsup1a.eps'))
saveas(gcf,fullfile(homepath,'newfigsup1a.png'))
savefig(fullfile(homepath,'newfigsup1a.fig'))

%save(fullfile(homepath,'fig1c'),'srs','ntrial','nsim','result','result2','rs','rs2','snr','snr2','M1s','M2s');
%% Figure2A
figure;
hold on;
idx=1:7; %[1,4,7]
idx_color=round(linspace(1,256,length(idx)));
color3=zeros(length(idx),3);
colorm=jet;
for i=1:length(idx)
    color3(i,:)=colorm(idx_color(i),:);
end

kk=1;
for i=idx
    X=(dWs{i}-0.5)/0.5;
    plot(diag(X),'Color',color3(kk,:),'LineWidth',1);
    %plot(diag(X),'LineWidth',1);
    kk=kk+1;
end
axis tight
colormap jet
colorbar
caxis([0.4 1.6])
xlabel('pre = post-synaptic firing rate (Hz)')
ylabel('change in synaptic strength')
box off
x0=100;
y0=200;
width=130.2;
height=102.8;

set(gcf,'units','points','position',[x0,y0,width,height])
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
saveas(gcf,fullfile(homepath,'fig2a.eps'))
saveas(gcf,fullfile(homepath,'fig2a.png'))
savefig(fullfile(homepath,'fig2a.fig'))
%% Figure2B
result=cell(length(dWs),1);
result2=cell(length(dWs),1);
result3=cell(length(dWs),1);
rs=cell(length(dWs),1);
rs2=cell(length(dWs),1);
rs3=cell(length(dWs),1);
snr=cell(length(dWs),1);
snr2=cell(length(dWs),1);
nsim=50;
ntrial=1;
srs=10:10:100;
M1s=cell(length(dWs),1);
M2s=cell(length(dWs),1);
for mm=1:8%idx
    X=dWs{mm,1};
    dWvar=dWvars{mm,1};
    [result{mm},result2{mm},result3{mm},rs{mm},rs2{mm},rs3{mm},snr{mm},snr2{mm},M1s{mm},M2s{mm}]=GPR_TE(X,dWvar,false,false,false,nsim,ntrial,srs);
end

save(fullfile(homepath,'fig2b'),'srs','ntrial','nsim','result','result2','rs','rs2','snr','snr2','M1s','M2s');
%%
figure;
tiledlayout(5,length(idx))
kk=1;
for i=idx
    a1=nexttile;
    imagesc('XData',f2,'YData',f2, 'CData',(dWs{i}-0.5)/0.5);
    hold on;
    plot(1:50,1:50,'--','Color',color3(kk,:));
    caxis([cl ch]) 
    colormap(a1,J)
    axis equal
    axis tight
    kk=kk+1;
end
colorbar
xlabel('presynaptic firing rate (Hz)')
ylabel('postsynaptic firing rate (Hz)')

for i=idx
    a2=nexttile;
    imagesc('XData',f2,'YData',f2, 'CData',(sqrt(dWvars{i})/0.5));
    caxis([0 0.1]) 
    colormap(a2,'gray')
    axis equal
    axis tight
end
colorbar
% xlabel('presynaptic firing rate (Hz)')
% ylabel('postsynaptic firing rate (Hz)')
    

for mm=idx
    mean_result=mean(squeeze(result{mm}),2);
    mean_result2=mean(squeeze(result2{mm}),2);
    mean_result3=mean(squeeze(result3{mm}),2);
    mean_snr=mean(squeeze(snr{mm}),2);
    mean_snr2=mean(squeeze(snr2{mm}),2);
    nexttile
    hold on;
%     l1=plot(srs,mean_result);
%     l2=plot(srs,mean_result2);
%     %l3=plot(srs,mean_result3);
%     plotmean(srs,squeeze(result{mm}),l1.Color);
%     plotmean(srs,squeeze(result2{mm}),l2.Color);
    %plotmean(srs,squeeze(result3{mm}),l3.Color);
    boxchart(squeeze(result{mm})')
    boxchart(squeeze(result2{mm})')
    xticklabels(10:10:100);
    axis tight
    %ylim([0.05,1])
    yline(0.2,'--')
end

xlabel('Number of observed points');
ylabel('Relative error');
legend('GPR','TE')

for mm=idx
    rtemp=squeeze(result2{mm});
    [~,medianindex]=min(abs(rtemp(2,:)-median(rtemp(2,:))));
    imgs=M2s{mm};
    a3=nexttile;
    imagesc('XData',f2,'YData',f2, 'CData',reshape(imgs{1,2,medianindex},[50,50]));
    caxis([cl ch]) 
    axis equal
    axis tight
    colormap(a3,J)
end

for mm=idx
    rtemp=squeeze(result2{mm});
    [~,medianindex]=min(abs(rtemp(10,:)-median(rtemp(10,:))));
    imgs=M2s{mm};
    a3=nexttile;
    imagesc('XData',f2,'YData',f2, 'CData',reshape(imgs{1,2,medianindex},[50,50]));
    caxis([cl ch]) 
    axis equal
    axis tight
    colormap(a3,J)
end

%% Figure 2 GPR vs TE Summary 
figure;
hold on;
for mm=1:7
    mean_result=mean(squeeze(result{mm}),2);
    mean_result2=mean(squeeze(result2{mm}),2);
    mean_snr=mean(squeeze(snr{mm}),2);
    mean_snr2=mean(squeeze(snr2{mm}),2);
    %scatter(result{mm}(:),result2{mm}(:),8,cc(:,mm)','filled');
    scatter(result{mm}(:),result2{mm}(:),8,snr2{mm}(:),'filled');
end
plot([0,1.2],[0,1.2],'--k')
xlabel('GPR relative error')
ylabel('TE relative error')
axis tight
caxis([0.1 0.8])
h=colorbar('XTick', 0.1:0.2:0.8);
colormap(gray*0.8)
ylabel(h, 'relative noise level')
x0=200;
y0=200;
width=130.2;
height=102.8;
set(gcf,'units','points','position',[x0,y0,width,height])
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
box off
saveas(gcf,fullfile(homepath,'fig2c.eps'))
saveas(gcf,fullfile(homepath,'fig2c.png'))
savefig(fullfile(homepath,'fig2c.fig'))
%% Figure 2C
idx=1:8;
resultc=cell(length(dWs),1);
resultc2=cell(length(dWs),1);
resultc3=cell(length(dWs),1);
rsc=cell(length(dWs),1);
rsc2=cell(length(dWs),1);
rsc3=cell(length(dWs),1);
snrc=cell(length(dWs),1);
snrc2=cell(length(dWs),1);
nsimc=50;
ntrialc=20;
srsc=20;
M1cs=cell(length(dWs),1);
M2cs=cell(length(dWs),1);
for mm=idx
    X=dWs{mm,1};
    dWvar=dWvars{mm,1};
    [resultc{mm},resultc2{mm},resultc3{mm},rsc{mm},rsc2{mm},rsc3{mm},snrc{mm},snrc2{mm},M1cs{mm},M2cs{mm}]=GPR_TE(X,dWvar,false,false,false,nsimc,ntrialc,srsc);
end
save(fullfile(homepath,'fig2c'),'srsc','ntrialc','nsimc','resultc','resultc2','rsc','rsc2','snrc','snrc2','M1cs','M2cs');
%%
load(fullfile(homepath,'fig2b'));
load(fullfile(homepath,'fig2c'));
%%
figure;
for mm=1:7
    nexttile
    hold on;
    mean_result=mean(squeeze(result{mm}),2);
    mean_result2=mean(squeeze(result2{mm}),2);
    mean_snr=mean(squeeze(snr{mm}),2);
    resultmatrix=squeeze(result{mm});
    
    mean_resultc=mean(squeeze(resultc{mm}),2);
    mean_resultc2=mean(squeeze(resultc2{mm}),2);
    mean_snrc=mean(squeeze(snrc{mm}),2);
    resultmatrixc=squeeze(resultc{mm});
    
    plot(20:20:100,mean_result(2:2:10),'--','Color',color3(mm,:));
    plot(20:20:100,mean_resultc(1:5),'Color',color3(mm,:));
    plotmean(20:20:100,resultmatrix(2:2:10,:),color3(mm,:));
    plotmean(20:20:100,resultmatrixc(1:5,:),color3(mm,:));
    axis tight;
end
xlabel('# of total observations');
ylabel('relative error');
legend('single Trials','repeated Trials')
saveas(gcf,fullfile(homepath,'figsup3e.eps'))
saveas(gcf,fullfile(homepath,'figsup3e.png'))
savefig(fullfile(homepath,'figsup3e.fig'))
%%
idx=1:7;
figure;
tile=tiledlayout(6,length(idx));
kk=1;
for i=idx
    a1=nexttile;
    imagesc('XData',f2,'YData',f2, 'CData',(dWs{i}-0.5)/0.5);
    hold on;
    plot(1:50,1:50,'--','Color',color3(i,:));
    caxis([cl ch]) 
    colormap(a1,J)
    axis equal
    axis tight
end
colorbar
xlabel('presynaptic firing rate (Hz)')
ylabel('postsynaptic firing rate (Hz)')

for i=idx
    a2=nexttile;
    imagesc('XData',f2,'YData',f2, 'CData',(sqrt(dWvars{i})/0.5));
    caxis([0 0.1]) 
    colormap(a2,'gray')
    axis equal
    axis tight
    if i~=idx(end)
        axis off
    end
end
colorbar
% xlabel('presynaptic firing rate (Hz)')
% ylabel('postsynaptic firing rate (Hz)')
    

for mm=idx
    mean_result=mean(squeeze(result{mm}),2);
    mean_result2=mean(squeeze(result2{mm}),2);
%    mean_result3=mean(squeeze(result3{mm}),2);
    mean_snr=mean(squeeze(snr{mm}),2);
    mean_snr2=mean(squeeze(snr2{mm}),2);
    nexttile
    hold on;
    l1=plot(srs,mean_result);
    l2=plot(srs,mean_result2);
%     %l3=plot(srs,mean_result3);
    plotmean(srs,squeeze(result{mm}),l1.Color);
    plotmean(srs,squeeze(result2{mm}),l2.Color);
    %plotmean(srs,squeeze(result3{mm}),l3.Color);
%     boxchart(squeeze(result{mm})')
%     boxchart(squeeze(result2{mm})')
%    xticklabels(10:10:100);
    plot(srs,mean_snr,'--k');
    axis tight
    %ylim([0.05,1])
    %yline(0.2,'--')
end

xlabel('number of observed points');
ylabel('relative error');
legend('GPR','TE')

% for mm=idx
%     rtemp=squeeze(result2{mm});
%     [~,medianindex]=min(abs(rtemp(2,:)-median(rtemp(2,:))));
%     imgs=M2s{mm};
%     a3=nexttile;
%     imagesc('XData',f2,'YData',f2, 'CData',reshape(imgs{1,2,medianindex},[50,50]));
%     %title(num2str(rtemp(2,medianindex)));
%     caxis([cl ch]) 
%     axis equal
%     axis tight
%     colormap(a3,J)
% 
% end

error1=zeros(length(idx),1);
for mm=idx
    rtemp=squeeze(result2{mm});
    [~,medianindex]=min(abs(rtemp(10,:)-median(rtemp(10,:))));
    imgs=M2s{mm};
    a3=nexttile;
    imagesc('XData',f2,'YData',f2, 'CData',reshape(imgs{1,10,medianindex},[50,50]));
    %title(num2str(rtemp(10,medianindex)));
    caxis([cl ch]) 
    axis equal
    axis tight
    colormap(a3,J)
    error1(mm)=rtemp(10,medianindex);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trials=1:20;
error2=zeros(length(idx),1);
for mm=idx
    rtemp=squeeze(resultc2{mm});
    [~,medianindex]=min(abs(rtemp(5,:)-median(rtemp(5,:))));
    imgs=M2cs{mm};
    ax4=nexttile;
    imagesc('XData',f2,'YData',f2, 'CData',reshape(imgs{5,1,medianindex},[50,50]));
    caxis([cl ch]) 
    %title(num2str(rtemp(5,medianindex)));
    axis equal
    axis tight
    colormap(ax4,J)
    error2(mm)=rtemp(5,medianindex);
    
end


for mm=idx
    mean_resultc=mean(squeeze(resultc{mm}),2);
    mean_resultc2=mean(squeeze(resultc2{mm}),2);
    %mean_result3=mean(squeeze(resultc3{mm}),2);
    mean_snrc=mean(squeeze(snrc{mm}),2);
    mean_snrc2=mean(squeeze(snrc2{mm}),2);
    nexttile
    hold on;
%     boxchart(squeeze(resultc{mm})')
%     boxchart(squeeze(resultc2{mm})')
%     yline(0.2,'--')

    l1=plot(trials,mean_resultc);
    l2=plot(trials,mean_resultc2);
    plotmean(trials,squeeze(resultc{mm}),l1.Color);
    plotmean(trials,squeeze(resultc2{mm}),l2.Color);
    plot(trials,mean_snrc,'--k');
    axis tight
    %ylim([0.05,0.9])
end
xlabel('number of trials');
ylabel('relative error');
legend('GPR','TE')

tile.TileSpacing = 'compact';
tile.Padding = 'compact';
%%
saveas(gcf,fullfile(homepath,'figsup2.eps'))
saveas(gcf,fullfile(homepath,'figsup2.png'))
savefig(fullfile(homepath,'figsup2.fig'))
%% fig 2C
figure;
hold on;
kk=1;
for mm=idx
    mean_result=mean(squeeze(resultc{mm}),2);
    mean_result2=mean(squeeze(resultc2{mm}),2);
    mean_result3=mean(squeeze(resultc3{mm}),2);
    mean_snr=mean(squeeze(snrc{mm}),2);
    mean_snr2=mean(squeeze(snrc2{mm}),2);
    scatter(mean_snr,mean_result,10,'filled','MarkerFaceColor',color3(kk,:));
    %plot(mean_snr,mean_result,'.','Color',color3(kk,:));
    %plotmean(mean_snr',squeeze(resultc{mm}),color3(kk,:));
    kk=kk+1;
end
x0=200;
y0=200;
width=130.2;
height=102.8;
set(gcf,'units','points','position',[x0,y0,width,height])
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
xlabel('relative noise level');
ylabel('relative error');
title('GPR');
saveas(gcf,fullfile(homepath,'fig2b.eps'))
saveas(gcf,fullfile(homepath,'fig2b.png'))
savefig(fullfile(homepath,'fig2b.fig'))

figure; 
hold on;
kk=1;
for mm=idx
    mean_result=mean(squeeze(resultc{mm}),2);
    mean_result2=mean(squeeze(resultc2{mm}),2);
    mean_result3=mean(squeeze(resultc3{mm}),2);
    mean_snr=mean(squeeze(snrc{mm}),2);
    mean_snr2=mean(squeeze(snrc2{mm}),2);
    scatter(mean_snr,mean_result2,10,'filled','MarkerFaceColor',color3(kk,:));
    kk=kk+1;
end
x0=200;
y0=200;
width=130.2;
height=102.8;
set(gcf,'units','points','position',[x0,y0,width,height])
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
xlabel('relative noise level');
ylabel('relative error');
title('TE');
saveas(gcf,fullfile(homepath,'figsup3a.eps'))
saveas(gcf,fullfile(homepath,'figsup3a.png'))
savefig(fullfile(homepath,'figsup3a.fig'))
%% fig 2C 2
figure;
hold on;
kk=1;
for mm=idx
    mean_result=mean(squeeze(resultc{mm}),2);
    mean_result2=mean(squeeze(resultc2{mm}),2);
    mean_result3=mean(squeeze(resultc3{mm}),2);
    mean_snr=mean(squeeze(snrc{mm}),2);
    mean_snr2=mean(squeeze(snrc2{mm}),2);
    scatter(squeeze(snrc{mm}),squeeze(resultc{mm}),10,'filled','MarkerFaceColor',color3(kk,:));
    kk=kk+1;
end
set(gcf,'units','points','position',[x0,y0,width,height])
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
xlabel('relative noise level');
ylabel('relative error');
axis tight
xlabel('relative noise level');
ylabel('relative error');
saveas(gcf,fullfile(homepath,'figsup3b.eps'))
saveas(gcf,fullfile(homepath,'figsup3b.png'))
savefig(fullfile(homepath,'figsup3b.fig'))

figure;
hold on;
kk=1;
for mm=idx
    mean_result=mean(squeeze(resultc{mm}),2);
    mean_result2=mean(squeeze(resultc2{mm}),2);
    mean_result3=mean(squeeze(resultc3{mm}),2);
    mean_snr=mean(squeeze(snrc{mm}),2);
    mean_snr2=mean(squeeze(snrc2{mm}),2);
    scatter(squeeze(snrc{mm}),squeeze(resultc2{mm}),10,'filled','MarkerFaceColor',color3(kk,:));
    kk=kk+1;
end
set(gcf,'units','points','position',[x0+100,y0,width,height])
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
xlabel('relative noise level');
ylabel('relative error');
axis tight
xlabel('relative noise level');
ylabel('relative error');
saveas(gcf,fullfile(homepath,'figsup3c.eps'))
saveas(gcf,fullfile(homepath,'figsup3c.png'))
savefig(fullfile(homepath,'figsup3c.fig'))
%%
figure;
hold on;
idx=1:7
for mm=idx
    mean_result=mean(squeeze(resultc{mm}),2);
    mean_result2=mean(squeeze(resultc2{mm}),2);
    mean_result3=mean(squeeze(resultc3{mm}),2);
    mean_snr=mean(squeeze(snrc{mm}),2);
    mean_snr2=mean(squeeze(snrc2{mm}),2);
    scatter(trials,mean_snr,10,'filled','MarkerFaceColor',color3(mm,:));
    temp(1)=mean_snr(1)
    for j=2:ntrialc
        temp(j)= temp(1)/sqrt(j);
    end
    plot(trials,temp,'Color',color3(mm,:))
end
axis tight;
xlabel('number of trials')
ylabel('noise level')
set(gcf,'units','points','position',[x0+100,y0,width,height])
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
saveas(gcf,fullfile(homepath,'figsup3d.eps'))
saveas(gcf,fullfile(homepath,'figsup3d.png'))
savefig(fullfile(homepath,'figsup3d.fig'))
%% Supp 1A: Rand vs Voronoi
result=cell(length(dWs),1);
result2=cell(length(dWs),1);
result3=cell(length(dWs),1);
rs=cell(length(dWs),1);
rs2=cell(length(dWs),1);
rs3=cell(length(dWs),1);
snr=cell(length(dWs),1);
snr2=cell(length(dWs),1);
nsim=50;
ntrial=1;
srs=10:10:100;
M1s=cell(length(dWs),1);
M2s=cell(length(dWs),1);
for mm=8
    X=dWs{mm,1};
    dWvar=dWvars{mm,1};
    [result{mm},result2{mm},result3{mm},rs{mm},rs2{mm},rs3{mm},snr{mm},snr2{mm},M1s{mm},M2s{mm}]=GPR_sampling(X,dWvar,false,false,false,nsim,ntrial,srs);
end
%% 
srs=10:10:100;
c1=[0 0 0];
c2=[0.5,0.5,0.5];
figure;
for mm=8
    mean_result=mean(squeeze(result{mm}),2);
    mean_result2=mean(squeeze(result2{mm}),2);
    mean_result3=mean(squeeze(result3{mm}),2);
    mean_snr=mean(squeeze(snr{mm}),2);
    mean_snr2=mean(squeeze(snr2{mm}),2);
    nexttile
    hold on;
    boxchart(squeeze(result{mm})','BoxFaceColor',c1,'MarkerColor',c1,'WhiskerLineColor',c1)
    boxchart(squeeze(result3{mm})','BoxFaceColor',c2,'MarkerColor',c2,'WhiskerLineColor',c2)
    
    %plotmean(srs,squeeze(result3{mm}),l3.Color);
    axis tight
end
xticklabels(srs);
plot(srs/10,mean_snr,'--k');
xlabel('number of observed points');
ylabel('relative error');
legend('Voronoi','Random')
x0=200;
y0=200;
width=130.2;
height=102.8;
%set(gcf,'units','points','position',[x0+100,y0,width,height])
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
saveas(gcf,fullfile(homepath,'figsup4a.eps'))
saveas(gcf,fullfile(homepath,'figsup4a.png'))
savefig(fullfile(homepath,'figsup4a.fig'))
%% Supp1B
result=cell(length(dWs),1);
result2=cell(length(dWs),1);
result3=cell(length(dWs),1);
rs=cell(length(dWs),1);
rs2=cell(length(dWs),1);
rs3=cell(length(dWs),1);
snr=cell(length(dWs),1);
snr2=cell(length(dWs),1);
nsim=50;
ntrial=20;
srs=20;
M1s=cell(length(dWs),1);
M2s=cell(length(dWs),1);
for mm=8
    X=dWs{mm,1};
    dWvar=dWvars{mm,1};
    [result{mm},result2{mm},result3{mm},rs{mm},rs2{mm},rs3{mm},snr{mm},snr2{mm},M1s{mm},M2s{mm}]=GPR_sampling(X,dWvar,false,false,false,nsim,ntrial,srs);
end
%%
trials=1:ntrial;

figure;
for mm=8
    mean_result=mean(squeeze(result{mm}),2);
    mean_result2=mean(squeeze(result2{mm}),2);
    mean_result3=mean(squeeze(result3{mm}),2);
    mean_snr=mean(squeeze(snr{mm}),2);
    mean_snr2=mean(squeeze(snr2{mm}),2);
    nexttile
    hold on;
%     l1=plot(trials,mean_result);
%     l2=plot(trials,mean_result2);
%    % l3=plot(trials,mean_result3);
%     plotmean(trials,squeeze(result{mm}),l1.Color);
%     plotmean(trials,squeeze(result2{mm}),l2.Color);
    boxchart(squeeze(result{mm})','BoxFaceColor',c1,'MarkerColor',c1,'WhiskerLineColor',c1)
    boxchart(squeeze(result3{mm})','BoxFaceColor',c2,'MarkerColor',c2,'WhiskerLineColor',c2)
    %plotmean(trials,squeeze(result3{mm}),l3.Color);
    axis tight
end
plot((1:1:20),mean_snr,'--k');
xlabel('number of trials');
ylabel('relative error');
legend('Voronoi','Random')
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
saveas(gcf,fullfile(homepath,'figsup4b.eps'))
saveas(gcf,fullfile(homepath,'figsup4b.png'))
savefig(fullfile(homepath,'figsup4b.fig'))
%%
lloydsAlgorithm(0.01*rand(20,1),zeros(20,1)+1/2, [0,0;0,1;1,1;1,0], 200, true)
title([])
xticks([0.2 0.4 0.6 0.8 1])
xticklabels(50*[0.2 0.4 0.6 0.8 1])
yticks([0.2 0.4 0.6 0.8 1])
yticklabels(50*[0.2 0.4 0.6 0.8 1])
set(gcf,'units','points','position',[x0+100,y0,width,height])
set(gca,'FontSize',8)
set(gca, 'FontName', 'Arial')
saveas(gcf,fullfile(homepath,'figsup4c.eps'))
saveas(gcf,fullfile(homepath,'figsup4c.png'))
savefig(fullfile(homepath,'figsup4c.fig'))

%% Rank of the 7 matrices (qualitatively or quantitatively)
rk_matrices=zeros(numel(idx),1);
rk_matrices2=zeros(numel(idx),1);
figure;
tiledlayout(numel(idx),3)
for mm=idx
    X=dWs{mm,1};
    ax1=nexttile;
    imagesc('XData',f2,'YData',f2, 'CData',((X-0.5)/0.5));
    axis equal; axis tight;
    caxis([cl ch]) 
    %title(num2str(rtemp(5,medianindex)));
    axis equal
    axis tight
    colormap(ax1,J)
    
    
    
    [u,s,v]=svd(X);
    diagval=diag(s);
    nl=mean(mean(squeeze(snr{mm}),2));
    rk=1;
    
    diagvall=zeros(size(diagval));
    diagvall(1:rk)=diagval(1:rk);
    sl=diag(diagvall);
    M_lr=u*sl*v';
    while norm(M_lr(:)-X(:))/norm(X(:))>0.01
        rk=rk+1;
        diagvall=zeros(size(diagval));
        diagvall(1:rk)=diagval(1:rk);
        sl=diag(diagvall);
        M_lr=u*sl*v';
    end
    rk_matrices(mm)=rk;
    ax2=nexttile;
     imagesc('XData',f2,'YData',f2, 'CData',((M_lr-0.5)/0.5));
    axis equal; axis tight;
    caxis([cl ch]) 
    %title(num2str(rtemp(5,medianindex)));
    axis equal
    axis tight
    colormap(ax2,J)
    title(['Rank = ', num2str(rk)])
    
    nexttile;
    plot(log(diag(s)));
    logs=log(diag(s));
    slope=zeros(49,1);
    for j=2:50
        slope(j)=logs(j)-logs(1);
    end
    %plot(slope);
    rk2=find(slope<-15,1);
    title(['Rank = ', num2str(rk2)])

    axis tight
    xlim([1 15])
end
box off
xlabel('rank')
ylabel('log(s)');
saveas(gcf,fullfile(homepath,'figsup5.eps'))
saveas(gcf,fullfile(homepath,'figsup5.png'))
savefig(fullfile(homepath,'figsup5.fig'))