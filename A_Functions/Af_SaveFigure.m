function Af_SaveFigure(figNumber,saveName,saveDir,varargin)
% Syntax: Af_SaveFigure(figNumber,saveName,saveDir,varargin)
%
% Function will save Figure figNumber in saveDir as a .fig and .png unless
% varargin specifies a different file type

if isempty(varargin)
    printString = '-dpng';
else
    printString = varargin{1};
end

figure(figNumber);

if strcmp(printString,'pdf')
    tightInset = get(gca, 'TightInset');
    position(1) = tightInset(1);
    position(2) = tightInset(2);
    position(3) = 1 - tightInset(1) - tightInset(3);
    position(4) = 1 - tightInset(2) - tightInset(4);
    set(gca, 'Position', position);
    saveas(gcf, [fullfile(saveDir,saveName),'.pdf']);
else
    
    % pos = get(gcf,'Position');
    % set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    savefig(fullfile(saveDir,saveName));
    set(gcf,'PaperPositionMode','auto');
    print(printString,fullfile(saveDir,saveName),'-r0');
end

