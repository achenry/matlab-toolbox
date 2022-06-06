
figure
Af_5MWSpdTqPlot

hold on
Power=2e5:2e5:6e6;

for p=1:length(Power)
    Ptq(p,:)=Power(p)./(W_hss(10:end)*pi/30);
    plot(W_hss(10:end),Ptq(p,:),'r');
end
axis([0,1400,0,5e4])
