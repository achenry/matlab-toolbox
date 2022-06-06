%function Shear = A_fMakeTVshear(TMax, DT,DC, freqHZ, order )
% Makes a time series shear, low pass filtered w/ order @ freq

% For testing:
TMax=600;
DT=.0125;
DC=.0;%.14;
freqHZ=.2;
order=2;
Amp=.2;
%XX Should I use rand or randn?
% Shear=randn(round(TMax/DT)+1,1);
% Shear=abs(DC*(Shear/2+1));
time=0:DT:TMax;

rs=.8*freqHZ*rand(1,5);

randSin=zeros(round(TMax/DT)+1,1);

N=5;
for n=1:N
    randSin=randSin+Amp/N*sin(2*pi*rs(n)*time'+rand);
end



Shear=2*Amp*(rand(round(TMax/DT)+1,1)-.5);

%Shear=abs(DC*(Shear/2+1));
%Shear=DC*(Shear/2+1);


Shear=Shear+randSin;
syms s
den=(1/(2*pi*freqHZ)*s+1)^order;
den=sym2poly(den);
[DFn,DFd]=tfdata(c2d(tf(1,den),DT));
DFn=cell2mat(DFn);
DFd=cell2mat(DFd);

Shear=filter(DFn,DFd,Shear);%,DC*10^-order*ones(1,order)');

maxi=max(Shear);
mini=min(Shear);

Shear=(Shear-mean(Shear))*2*Amp/(maxi-mini)+DC;

save('Ashear.mat','Shear');
figure
plot(time',Shear)
%end

