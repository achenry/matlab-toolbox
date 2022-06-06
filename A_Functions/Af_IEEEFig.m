function h = Af_IEEEFig(h,width,varargin)
% Varargin(1) = FontSize
% varargin(2) = figure height

if numel(varargin) > 0
    figureFontSize = varargin{1};
else
    figureFontSize = [];
end

if isempty(figureFontSize)
    figureFontSize = 9;  %for wind energy science
else
    figureFontSize = varargin{1};
end
% 
set(h,'defaulttextinterpreter','latex');
set(h, 'defaultAxesTickLabelInterpreter','latex');
set(h, 'defaultLegendInterpreter','latex');

% set(0,'defaulttextinterpreter','none');
% set(0, 'defaultAxesTickLabelInterpreter','none');
% set(0, 'defaultLegendInterpreter','none');

% set(groot, 'defaultAxesFontName','Times New Roman')
% set(groot, 'defaultTextFontName','Times New Roman');

%% Set Font Sizes

ha = get(gcf,'Children');
for ih = 1:length(ha)
    set(ha(ih), 'FontSize', figureFontSize);
    set(ha(ih), 'FontName', 'Times New Roman');
    
    if strcmp(get(ha(ih),'type'),'axes')
        set(ha(ih), 'TitleFontSizeMultiplier', 1);
    end
end

%% Size Image

set(gcf,'Units','inches');
pos = get(h,'Position');
baseWidth = 3;

width = width * baseWidth;
if numel(varargin) > 1
    height = varargin{2};
else
    height = pos(4);
end

fig = gcf;
fig.Units = 'inches';
fig.Position = [5,0,width,height];
fig.PaperSize = [width,height];


% set(gcf,...
%     'Units','inches',...
%     'PaperSize',[width,height],...
%     'PaperPosition',[0,0,width,height],...
%     'Position',[5,1,width,height]);