function plotmean(x,MSE,color)
    y=mean(MSE,2);
    e=std(MSE,0,2)./sqrt(size(MSE,2));
    yu=y+e;
    yl=y-e;
%     plot(x,y,color);
%     hold on;
    area=fill([x fliplr(x)], [yu' fliplr(yl')], color, 'linestyle', 'none');
    area.FaceAlpha=0.2;
end