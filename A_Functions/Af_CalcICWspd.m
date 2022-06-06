function [ IC_WSPD ] = Af_CalcICWspd( IC_Wr,IC_Bp,IC_Tg,IC_Wrdot)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global WT_DATA WT_XDATA WT_YDATA TSRtest rho Rrad Ngear Jtot

%0=0.5*rho*pi*Rrad^2*V^3*CP-Ngear*X(3)
CPslice=interp2(WT_XDATA,WT_YDATA,WT_DATA{2},IC_Bp,TSRtest);
%plot(TSRtest,CPslice)
WStest=IC_Wr*Rrad./TSRtest;

EQ=0.5*rho*pi*Rrad^2*WStest.^3.*CPslice'-Ngear*IC_Tg-IC_Wrdot*Jtot;
%figure, plot(WStest,EQ);

IC_WSPD=interp1(EQ,WStest,0,'linear','extrap');

end

