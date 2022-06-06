function A_fMakeSinWind(dcwind, freq, amp, Tmax, dt, shear)
%% Make Wind File
% Jacob Aho
% EXAMPLE
% ! Steady wind file created 7/22/98 7:58:37 PM by YawDynVB Version 2.0
% ! Time	Wind	Wind	Vert.	Horiz.	Vert.	LinV	Gust
% !	        Speed	Dir	    Speed	Shear	Shear	Shear	Speed
%   0.0	     30	    0       -1      0.1     0.14	0       0
 
samples=round(Tmax/dt+1);

W=zeros(samples, 8);

W(:,1)=0:dt:Tmax;

W(:,2)=dcwind+amp*sin(2*pi*W(:,1)*freq);

%W(samples/2:samples*3/4,2)=18;

W(:,4)=0;
W(:,5)=0;
W(:,6)=shear;
W(:,7)=0;
W(:,8)=0;

dlmwrite('WindData\sine.wnd', W, ' ')
disp(['Sinusoidal Wind File Made: DC=', num2str(dcwind), ' m/s, Freq=', num2str(freq), ' Amp=', num2str(amp), ' Shear=', num2str(shear)])