function [Wtype,varargout]=Af_MakeWind(WindDir,Disturbance,Simulation,Wind_Params,Mplot,PARSIM)
% function [Wtype,varargout]=Af_MakeWind(WindDir,Disturbance,Simulation,TYPE,TMax, DT, WD, VS, Hshear, Vshear, LVshear, GS, Wind_Params,Mplot)
% J.Aho 9/22/11
% This function will make a 'User Defined' wind file for use with FAST
% The wind file is generated based on the input to the function
% TYPE:
%
% 9/6/17: Added Peak load wind type
% 12/4/17: streamlined inputs

TYPE    = Disturbance.Type;
TMax    = Simulation.TMax;
DT      = Simulation.DT;

if TYPE==1 % Constant Wind
    %WindParams(1)= [DCvalue]
    Wtype='const';
    
elseif TYPE==2  % Stepped input from Vmin to Vmax then back down in steps of dWind
    Wtype='step';
    %WindParams(1:3)=[Vmin, Vmax, dWind]
    
elseif TYPE==3  % Ramp wind going from Vmin to Vmax to min, ramp starts at Stime
    Wtype='ramp';
    %WindParams2(1:4)=[min, max, Ramp Start Time, Ramp Duration]
    
elseif TYPE==4   % Sinusoidal wind of DC+amp*sin(2*pi*freq*t)
    Wtype='sine';
    %WindParams2(1:4)=% [DC, freq (Hz), amp]
elseif TYPE==5
    Wtype='sstep';
elseif TYPE==6
    Wtype='yaw_step';
    
elseif TYPE==7
    Wtype='EOG';
    
elseif TYPE==8
    Wtype='ECD';
    
elseif TYPE == 9
    Wtype = 'PeakTest';
    
elseif TYPE == 10
    Wtype = 'ConstRake';
    
end

% The following parameters are fully described in the TurbSim/Aerodyn docs.

% TMax= Max time, DT=Time step, WD= Wind Direction, VS= Vertical speed
% Hshear=Linear horizontal shear, Vshear= Power Law Vertical Shear,
% LVshear= Linear Vertical Shear, GS= Gust Speed
% Mplot==1 will produce a plot of the generated wind file (hub height)

% Note: If any shear parameter input is a scalar, it will produce a wind
% file with a constant shear.  If the input shear parameter is a vector,
% it must be the same length as the generated wind file and is used as the
% time varying shear.

% EXAMPLE with Header
% ! Time	Wind	Wind	Vert.	Horiz.	Vert.	LinV	Gust
% !	        Speed	Dir	    Speed	Shear	Shear	Shear	Speed
%   0.0	     30	    0       -1      0.1     0.14	0       0

% Generated wind file is stored in the PWD\WindData folder with
% the appropriate name for the type.



p1=Wind_Params(1);
p2=Wind_Params(2);
p3=Wind_Params(3);

if length(Wind_Params == 4)
    p4=Wind_Params(4);
    
elseif length(Wind_Params == 5)
    p4=Wind_Params(4);
    p5=Wind_Params(5);
    
elseif length(Wind_Params == 6)
    p4=Wind_Params(4);
    p5=Wind_Params(5);
    p6=Wind_Params(6);
end

samples=round(TMax/DT+1);
W=zeros(samples, 8);

W(:,1)=0:DT:TMax;

WD = 0;
if length(WD)>1
    W(:,3)=WD;
else
    W(:,3)=ones(samples,1)*WD;
end

VS = 0; %no vertical speed for now
if length(VS)>1
    W(:,4)=VS;
else
    W(:,4)=ones(samples,1)*VS;
end

Hshear = Disturbance.LHshear;
if length(Hshear)>1
    W(:,5)=Hshear;
else
    W(:,5)=ones(samples,1)*Hshear;
end

Vshear = Disturbance.Vshear;
if length(Vshear)>1
    W(:,6)=Vshear;
else
    W(:,6)=ones(samples,1)*Vshear;
end

LVshear = Disturbance.LVshear;
if length(LVshear)>1
    W(:,7)=LVshear;
else
    W(:,7)=ones(samples,1)*LVshear;
end

GS = 0; %gust = 0 for now
if length(GS)>1
    W(:,8)=GS;
else
    W(:,8)=ones(samples,1)*GS;
end


switch Wtype
    case 'const'  % p1=DC VALUE
        W(:,2)=ones(samples,1)*p1;
    case 'step' %p1=min, p2=max, p3=ds
        
        p2=p2+p3;
        
        sgn=sign(p2-p1);
        
        W(:,2)=ones(samples,1)*p1;
        steptime=round(samples/2/(sgn*(p2-p1)/p3));
        
        for n=1:2*sgn*(p2-p1)/p3
            b=(n-1)*steptime+1;
            e=n*steptime+1;
            if e>samples;
                e=samples;
            end
            if n<=sgn*(p2-p1)/p3
                W(b:e,2)=p1+sgn*(n-1)*p3;
            end
            if n>sgn*(p2-p1)/p3
                W(b:e,2)=p2-sgn*(n-sgn*(p2-p1)/p3)*p3;
            end
        end
        
        
    case 'sine' % p1=dcwind, p2=freq, p3=amp,
        
        
        W(:,2)=p1+p3*sin(2*pi*W(:,1)*p2);
        
    case 'ramp' % p1=min, p2=max, p3=ST
        
        n=1;
        W(1:samples,2)=p1;
        W(round(p3/DT)+1:round((p3+p4/2)/DT)+1,2)=(p2-p1)/((p4)/2).*(W(round(p3/DT)+1:round((p3+p4/2)/DT+1),1)-p3)+p1;
        W(round((p3+p4/2)/DT)+1:round((p3+p4)/DT)+1,2)=(p1-p2)/((p4)/2).*(W(round(((p3+p4+p3)/DT)/2)+1:round((p3+p4)/DT)+1,1)-(p4+p3+p3)/2)+p2;
        
    case 'sstep' % p1:4  = [ Init. Speed   Ramp Start      Final Speed     Ramp End   ]
        tt = W(:,1);
        if p2 == p4     %just a step
            vv = p1*ones(samples,1);
            vv(tt>p2) = p3;
        else
            X = [0,p2,p4,TMax];
            V = [p1,p1,p3,p3];
            vv = interp1(X,V,tt);
        end
        W(:,2) = vv;
        
    case 'yaw_step' % p1:5  [ Wind Speed    Initial Yaw     Start Time      Final Yaw       End Time ]
        
        tt = W(:,1);
        
        if p3 == p5     %just a step
            yy = p2*ones(samples,1);
            yy(tt>p3) = p4;
            
        else            %ramp
            X = [0,p3,p5,TMax];
            V = [p2,p2,p4,p4];
            yy = interp1(X,V,tt);
        end
        W(:,2) = ones(samples,1)*p1./cosd(yy)+1.373157942907710e-08;
        W(:,2) = ones(samples,1)*p1;
        W(:,3) = yy;
        
        
    case 'EOG'
        tt = W(:,1);
        tt_gust = (0:DT:p3)';
        U  = p1;
        gust = [zeros(p2/DT,1);
            -0.37*p4*sin(3*pi*tt_gust/p3).*(1-cos(2*pi*tt_gust/p3));
            zeros(samples-p2/DT-p3/DT-1,1)];
        W(:,2) = U+gust;
        
    case 'ECD'
        n=1;
        W(1:samples,2)=p1;
        W((p2/DT)+1:((p2+p4)/DT)+1,2)= p1+0.5*15*(1-cos(pi*[((0)/DT)+1:((p4)/DT)+1]/(p4/DT)));
        W(round((p2+p4)/DT)+1:end,2)= p3;
        if p1<4
            theta_cg=180;
        else
            theta_cg=720/p1;
        end
        W(1:samples,3) = 0;
        W((p2/DT)+1:((p2+p4)/DT)+1,3) = 0.5*theta_cg*(1-cos(pi*[((0)/DT)+1:((p4)/DT)+1]/(p4/DT)));
        W(round((p2+p4)/DT)+1:end,3)= theta_cg;
        
        
    case 'PeakTest'
        startTime = 200;                    %start time of first gust
        extraTime = 100;                    %time after test
        periods = [60,45,30,20,10,5];       %sinusoid periods
        
        DT = 1/80;
        T = periods(1);
        tt = 0:DT:TMax;
        gust = zeros(length(periods),length(tt));
        
        startGust = [startTime,startTime + cumsum(periods(1:end-1))];
        
        Amp = 6;                            %gust amplitude
        startU = 9.3;                       %wind speed to start gust
        
        for iGust = 1:size(gust,1)-1
            T = periods(iGust);
            gust(iGust,:) = (Amp/2 * (-cos(2*pi/T * (tt - startGust(iGust)))+1)).*(tt >= startGust(iGust) & tt <= startGust(iGust+1));
        end
        
        wind = sum(gust,1) + startU;
        W(:,1) = tt;
        W(:,2) = wind;
        
        if 1
            figure(100);
            plot(tt,wind);
        end
        
    case 'ConstRake' % p1=min, p2=t_start, p3=t_end, p4 = du
        
        tt      = 0:DT:TMax;
        b       = p1 - p4/20 * p2;
        m       = p4 /20;
        maxWS   = 25;
        tend    = (maxWS-b)*(20/p4);
        
        W(:,2)      = p1;
        W(tt>p2,2)  = m * tt(tt>p2) + b;
        W(tt>tend,2)    = maxWS;
        
        if 1
            figure(100);
            plot(tt,W(:,2));
        end
end

varargout{1} = W;
dlmwrite(fullfile(WindDir,[Wtype,PARSIM,'.wnd']), W, 'delimiter',' ','precision',8)
disp([Wtype,' Wind File Made'])
if Mplot==1
    figure
    plot(W(:,1),W(:,2))
    xlabel('Time [s]')
    ylabel('HH Velocity [m/s]')
    grid on
end