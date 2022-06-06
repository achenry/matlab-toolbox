function Af_PlotFig(linewidth,fontsize)

% Fontsize as a vector: Font, XYZlabels, Title

%% For Figures

axesHandles = get(gcf,'Children');

for iAxes = 1:length(axesHandles)

if fontsize>0
    if size(fontsize)==1
    set(axesHandles(iAxes),'FontSize',fontsize)
    set(get(axesHandles(iAxes),'XLabel'),'FontSize',fontsize)
    set(get(axesHandles(iAxes),'YLabel'),'FontSize',fontsize)
    set(get(axesHandles(iAxes),'ZLabel'),'FontSize',fontsize)
    set(get(axesHandles(iAxes),'Title'),'FontSize',fontsize)
    else
       set(axesHandles(iAxes),'FontSize',fontsize(1))
    set(get(axesHandles(iAxes),'XLabel'),'FontSize',fontsize(2))
    set(get(axesHandles(iAxes),'YLabel'),'FontSize',fontsize(2))
    set(get(axesHandles(iAxes),'ZLabel'),'FontSize',fontsize(2))
    set(get(axesHandles(iAxes),'Title'),'FontSize',fontsize(3))

    end
end
if linewidth>0
    hls123 = findobj(gcf,'Type','Line');
    set(hls123,'LineWidth',linewidth)
end
end