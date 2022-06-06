%add in pitch model and do mbc:

%get system input and output names
avg_Bdsys=Fgui_data.AVG_BDSYS;
sys_outputs=avg_Bdsys.outputnames;
sys_inputs=avg_Bdsys.inputnames;

%get models at all rotor positions
A=Fgui_data.AMat;
B=Fgui_data.BMat;
Bd=Fgui_data.BdMat;
C=Fgui_data.CMat;
D=Fgui_data.DMat;
Dd=Fgui_data.DdMat;

%make sure system has inputs in desired order
Bdsys_azs=ss(A,[B,Bd],C,[D,Dd]);
Bdsys_azs.inputnames=sys_inputs;
% % Bdsys_azs=Bdsys_azs(:,{'P1','P2','P3','HHwind','VLshear','Hshear'});
sys_inputs=Bdsys_azs.inputnames;

%get system size number of rotor positions
[p,m,n_theta]=size(Bdsys_azs);
n=order(Bdsys_azs);n=n(1);
n=n+3;      %will add three actuator states
theta=linspace(0,2*pi/n_theta*(n_theta-1),n_theta);
omega=Fgui_data.NominalSpeed/60*2*pi;

%initialize matrices for avg models
p=p+6; %will add pitch and pitch rate outputs
aMBC=zeros(n,n);aAvg=aMBC;
bMBC=zeros(n,m);bAvg=bMBC;
cMBC=zeros(p,n);cAvg=cMBC;
dMBC=zeros(p,m);dAvg=dMBC;

%make actuator model with pitch and pitch rate outputs
if ~exist('actpole')
    actpole=30;
end;
actmod=ss(-actpole,actpole,[1;-actpole],[0;actpole]);
actmod=append(actmod,actmod,actmod);

Nout=length(sys_outputs);
Nin=length(sys_inputs);
new_inputs={sys_inputs{:},'pr1','pr2','pr3'};
new_outputs={sys_outputs{:},'P1','P2','P3','pr1','pr2','pr3'};
MBC2SHEAR=diag([1,mbc2Vshear,mbc2Hshear]);

for index=1:n_theta
    Bdsys=Bdsys_azs(:,:,index);
    
    %add feedthrough to pass pitchrate and pitch from actuator model
    Bdsys=append(Bdsys,1,1,1);
    Bdsys=Bdsys*[eye(Nin);[zeros(3,1),eye(3),zeros(3,Nin-4)]]; %feed pitch to output
    Bdsys=append(Bdsys,1,1,1);
    Bdsys.inputnames=new_inputs;
    Bdsys.outputnames=new_outputs;

    %splice pitch and pitch rate inputs
    Bdsys=Bdsys(:,{'Tgen','P1','pr1','P2','pr2','P3','pr3','HHwind','VLshear','Hshear'});
    
    %insert pitch actuator models
    turbine=Bdsys*append(1,actmod,1,1,1);
    turbine.inputnames=sys_inputs;
    
    %form MBC atoms for use in transformation
    MBCa=diag([1/3,2/3,2/3])*...
            [ 1                 1                         1
              cos(theta(index)) cos(theta(index)+2*pi/3)  cos(theta(index)+4*pi/3)
              sin(theta(index)) sin(theta(index)+2*pi/3)  sin(theta(index)+4*pi/3) ];
    %inverse 
    MBCai=  [ 1                 cos(theta(index))         sin(theta(index)) 
              1                 cos(theta(index)+2*pi/3)  sin(theta(index)+2*pi/3)
              1                 cos(theta(index)+4*pi/3)  sin(theta(index)+4*pi/3) ];
    %d/dt inverse
    MBCaip= [ 0                 -sin(theta(index))        cos(theta(index)) 
              0                 -sin(theta(index)+2*pi/3) cos(theta(index)+2*pi/3)
              0                 -sin(theta(index)+4*pi/3) cos(theta(index)+4*pi/3) ]*omega;
    %d^2/dt^2 inverse
    MBCaipp=[ 0                 -cos(theta(index))        -sin(theta(index)) 
              0                 -cos(theta(index)+2*pi/3) -sin(theta(index)+2*pi/3)
              0                 -cos(theta(index)+4*pi/3) -sin(theta(index)+4*pi/3) ]*omega^2;

    %build transformations for state displacements:
    MBCistate=blkdiag(1,1,MBCai);               %T:  DISP_nr => DISP
    MBCstate=blkdiag(1,1,MBCa);                 %Ti: DISP    => DISP_nr
    
    %build partial transformation for state velocities:
    MBCistatep=blkdiag(0,0,MBCaip);             %dT/dtheta
    MBCistatepp=blkdiag(0,0,MBCaipp);           %d^2T/dtheta^2
    MBCi1=[ MBCistate   zeros(5,5) zeros(5,3)
            MBCistatep  MBCistate  zeros(5,3)
            zeros(3,5)  zeros(3,5) MBCai     ];
        
    MBCi2=[ MBCistatep  zeros(5,5)   zeros(5,3)
            MBCistatepp 2*MBCistatep zeros(5,3)
            zeros(3,5)  zeros(3,5)   MBCaip     ];
        
% %     MBCi1=[ MBCistate   zeros(5,5) zeros(5,3) zeros(5,3)
% %             MBCistatep  MBCistate  zeros(5,3) zeros(5,3)
% %             zeros(3,5)  zeros(3,5) MBCai      zeros(3,3)
% %             zeros(3,5)  zeros(3,5) zeros(3,3) MBCai     ];
% %         
% %     MBCi2=[ MBCistatep  zeros(5,5)   zeros(5,3) zeros(5,3)
% %             MBCistatepp 2*MBCistatep zeros(5,3) zeros(5,3)
% %             zeros(3,5)  zeros(3,5)   MBCaip     zeros(3,3)
% %             zeros(3,5)  zeros(3,5)   zeros(3,3) MBCaip    ];

    Tbari=blkdiag(MBCstate,MBCstate,MBCa);
    
    %Input and output transformations
    MBCinput=blkdiag(1,MBCai,eye(3));         %TI:  NR     => non-NR
    MBCoutput=blkdiag(1,MBCa,MBCa,MBCa);      %TOi: non-NR => NR
    
    %apply MBC transformation
    [a,b,c,d]=ssdata(turbine);
    ambc=Tbari*(a*MBCi1-MBCi2);
    bmbc=Tbari*b*MBCinput;
    cmbc=MBCoutput*c*MBCi1;
    dmbc=MBCoutput*d*MBCinput;
    
    %accumulate state matrices to form avg model
    %TI:  blade local  => NR shear
    BLOCALinput=blkdiag(eye(4),MBC2SHEAR*MBCa);  
    bblocal=b*BLOCALinput; %shear in = MBC2SHEAR*MBCa * [w1,w2,w3]'
    dblocal=d*BLOCALinput;
    aAvg=aAvg+a;
    bAvg=bAvg+bblocal;
    cAvg=cAvg+c;
    dAvg=dAvg+dblocal;
    aMBC=aMBC+ambc;
    bMBC=bMBC+bmbc;
    cMBC=cMBC+cmbc;
    dMBC=dMBC+dmbc;
end;
aAvg=aAvg/n_theta;
bAvg=bAvg/n_theta;
cAvg=cAvg/n_theta;
dAvg=dAvg/n_theta;
aMBC=aMBC/n_theta;
bMBC=bMBC*blkdiag(eye(4),MBC2SHEAR)/n_theta; %shear in = MBC2SHEAR * [wC,wV,wH]'
cMBC=cMBC/n_theta;
dMBC=dMBC*blkdiag(eye(4),MBC2SHEAR)/n_theta; %shear in = MBC2SHEAR * [wC,wV,wH]'
turbine=ss(aAvg,bAvg,cAvg,dAvg);

%drop rotor position state
turbine=modred(turbine,1,'truncate');
turbine.StateNames={ 'dtrain','flap1','flap2','flap3','gen',...
                            'dtraindt','flap1v','flap2v','flap3v',...
                            'act1','act2','act3'};
turbine.outputnames={sys_outputs{:},'P1','P2','P3','pr1','pr2','pr3'};
turbine.inputnames={sys_inputs{1:4},'w1','w2','w3'};

turbineMBC=ss(aMBC,bMBC,cMBC,dMBC);
turbineMBC.outputnames={'HSShftV','RootMyC','RootMyV','RootMyH','pC','pV','pH','prC','prV','prH'};
turbineMBC.inputnames={'Tgen','pC','pV','pH','wC','wV','wH'};
turbineMBC=modred(turbineMBC,1,'truncate'); %drop rotor position state
turbineMBC.StateNames={ 'dtrain','flapC','flapV','flapH','gen',...
                            'dtraindt','flapdtC','flapdtV','flapdtH',...
                            'actC','actV','actH'};

