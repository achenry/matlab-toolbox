function Af_PlotFig_Pres(linewidth,fontsize)



%% For Figures
if fontsize>0
    set(gca,'FontSize',fontsize)
    set(get(gca,'XLabel'),'FontSize',fontsize)
    set(get(gca,'YLabel'),'FontSize',fontsize)
    set(get(gca,'ZLabel'),'FontSize',fontsize)
    set(get(gca,'Title'),'FontSize',fontsize)
end
if linewidth>0
    hls123 = findobj(gcf,'Type','Line');
    set(hls123,'LineWidth',linewidth)
end
end