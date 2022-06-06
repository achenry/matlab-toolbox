% Gen Speed vs Tq plot

hold on

W_hss=0:.1:(1.2*WratedRPM*Ngear);

plot([WratedRPM*Ngear,WratedRPM*Ngear],[0,1.2*Trated],'r--')
hold on
plot([0,1.2*WratedRPM*Ngear],[Trated,Trated],'r--')
plot(W_hss,K*W_hss.^2,'g','LineWidth',3)


title('Gen Torque Vs. Speed')
xlabel('Gen Speed (RPM)')
ylabel('Gen Torque (Nm)')
xlims=[-inf inf];
ylims=[-inf, inf];
axis([xlims ylims])
Af_PlotFig(0,14);
grid on
