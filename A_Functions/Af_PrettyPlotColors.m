function Af_PrettyPlotColors(p,varargin)
%set color...fancy!
% 
% Input: p  plot handles to set color of

if ~isempty(varargin)
    colorMap    = colormap(varargin{1});
else
    colorMap    = colormap('jet');
end
colorInterval   = floor(size(colorMap,1)/length(p));        %evenly distribute

colorIndex      = 1:colorInterval:length(p)*colorInterval+1;
for iPlot = 1:length(p)
    p(iPlot).Color = colorMap(colorIndex(iPlot),:);
end