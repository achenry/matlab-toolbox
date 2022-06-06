%% Plot Platform Freqs

xlim([0.01,2]);

yl = get(gca,'YLim');

options = bodeoptions;
options.FreqUnits = 'Hz';

freqs = [0.036,0.502];

hold on;
loglog(repmat(freqs,[2,1]),repmat(yl,[length(freqs),1])','k')
% loglog(PP.ATLAS.Platform.Freq,PP.ATLAS.Platform.SpecA,'x','LineWidth',2,'Color',h.Color,'MarkerSize',10);
hold off;

grid on;