%% Progress Bar Script
% not a function, need nSims, iSim

% progress bar
if ~exist('hProgBar')
    hProgBar = figure(1000);
elseif ~isvalid(hProgBar)
    hProgBar = figure(1000);
else
    set(0,'CurrentFigure',hProgBar);
end
bar(iSim/nSims);
ylim([0,1]);
drawnow;