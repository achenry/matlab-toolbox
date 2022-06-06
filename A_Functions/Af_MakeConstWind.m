function  A_fMakeConstWind(c,  Tmax, dt, shear)
%% Make a constant wind file wind file
% Jacob Aho

% EXAMPLE of what columns are!
% ! Steady wind file created 7/22/98 7:58:37 PM by YawDynVB Version 2.0
% ! Time	Wind	Wind	Vert.	Horiz.	Vert.	LinV	Gust
% !	        Speed	Dir	    Speed	Shear	Shear	Shear	Speed
%   0.0	     30	    0       -1      0.1     0.14	0       0
 

samples=round(Tmax/dt+1);

W=zeros(3, 8);

W(:,1)=[0,dt,Tmax]';

W(:,2)=[c c c]';

W(:,4)=0;
W(:,5)=0;
W(:,6)=[shear, shear, shear]';
W(:,7)=0;
W(:,8)=0;


dlmwrite('FAST_IF\WindData\const.wnd', W, ' ')

disp(['Constant Wind File Made at ', num2str(c), ' m/s, shear=', num2str(shear)])