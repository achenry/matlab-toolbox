function [hpL,hpB] =  Af_plotMinPitch(UU,PP,HOLD,figNum,lineStyle)

if nargin < 4
    figNum = 1050;
end

if nargin < 5
    lineStyle = '-';
end

figure(figNum);
uu = linspace(5,25);

pp = interp1(UU,PP,uu,'pchip');

if HOLD, hold on; end
hpL = plot(uu,pp,'LineWidth',2,'LineStyle',lineStyle);

ind = get(gca,'ColorOrderIndex');
if ind == 1
    set(gca,'ColorOrderIndex',7);
else
    set(gca,'ColorOrderIndex',ind-1);
end

hold on;
hpB = plot(UU,PP,'x','LineWidth',2,'MarkerSize',10);
hold off;
xlim([5,25]);
grid on;

xlabel('Wind Speed (m/s)');
ylabel('Min. Blade Pitch (deg.)');