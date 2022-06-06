% Gen Speed vs Tq plot

hold on


W_hss=0:.1:1400;


% Define Constants
m1_5 = 96.5338;         %Region 1.5 slope (Nm/rpm)
c1_5 = -64677.65123;    %Region 1.5 offset (Nm)
K = 0.025576386;        %Region 2 gain (Nm/rpm^2)
m2_5= 412.076;          %Region 2.5 slope (Nm/rpm)
c2_5 = -435288.3165;     %Region 2.5 offset (Nm)
P0 = 5e6/0.944;         %Region 3 raw power (W)

for count=1:14001
if (W_hss(count) >= 1161.9632) %Region 3 calcs
        genTorque(count) = 30/pi * P0 / W_hss(count);
elseif (W_hss(count) >= 1136.4978)                %Region 2.5
    genTorque(count) = m2_5 * W_hss(count) + c2_5;
elseif (W_hss(count) >= 871)                    %Region 2
    genTorque(count) = K * W_hss(count)^2;
elseif (W_hss(count) > 670)                    %Region 1.5
    genTorque(count) = m1_5 * W_hss(count) + c1_5;
else                                    %Region 1
    genTorque(count) = 0;
end
end

plot(W_hss,genTorque,'LineWidth',2)
hold on
plot(W_hss,K * W_hss.^2,'K','LineWidth',2)
title('Gen Spd Vs. Torque- Baseline')
xlabel('Gen Speed (RPM)')
ylabel('Gen Torque (Nm)')
xlims=[-inf inf];
ylims=[-inf, inf];
axis([xlims ylims])
grid on
