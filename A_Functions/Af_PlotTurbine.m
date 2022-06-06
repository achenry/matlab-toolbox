function Af_PlotTurbine(az,yaw,varargin)
% Af_PlotTurbine(azimuith,yaw) will visualize the NREL 5MW turbine at a
% specified azimuith and yaw angle in degrees

%%%%% TODO %%%%%%
% 1. Rotor Coning and Tilt


%Load 5MW Turbine Params
C_5MW_Params;

%% Tower

Tower = [0,0;
    0,0;
    0,87];

figure(90);
if ~isempty(varargin)
    eval(varargin{1});
end
plot3(Tower(1,:),Tower(2,:),Tower(3,:),'LineWidth',7); hold on;
grid on;
axis([-100,100,-100,100,0,200])
zlabel('z');xlabel('x');ylabel('y');

%% Blades


R = Parameters.Turbine.Rrad;

Blade1 = [Tower(:,2),Tower(:,2) + R*[sind(az)*cosd(yaw) ; sind(az)*sind(yaw) ; cosd(az)]];
Blade2 = [Tower(:,2),Tower(:,2) + R*[sind(az+120)*cosd(yaw) ; sind(az+120)*sind(yaw) ; cosd(az+120)]];
Blade3 = [Tower(:,2),Tower(:,2) + R*[sind(az-120)*cosd(yaw) ; sind(az-120)*sind(yaw) ; cosd(az-120)]];


plot3(Blade1(1,:),Blade1(2,:),Blade1(3,:),'LineWidth',3); hold on;
plot3(Blade2(1,:),Blade2(2,:),Blade2(3,:),'LineWidth',3); hold on;
plot3(Blade3(1,:),Blade3(2,:),Blade3(3,:),'LineWidth',3); hold on;


%% Rotor Disk

Rotor = repmat(Tower(:,2),1,100) +  R*[sind(linspace(0,360))*cosd(yaw);
    sind(linspace(0,360))*sind(yaw);
    cosd(linspace(0,360))];

rot = fill3(Rotor(1,:),Rotor(2,:),Rotor(3,:),[.9,.9,.9]);
set(rot,'FaceAlpha',.4);



%% Wind Vector

if 0
    r = .75*R;
    
    section = Blade1(:,1) + r*[sind(az)*cosd(yaw) ; sind(az)*sind(yaw) ; cosd(az)];
    
    U = 40;
    
    U_inf = [section - [0;U;0] , section];
    
    th_ma1 = yaw*cosd(az);
    th_mb1 = yaw*sind(az);
    
    u1_a = U*cosd(th_ma1)*cosd(th_mb1);
    u1_t = -U*sind(th_ma1)*cosd(th_mb1);
    
    U1_a = [section - u1_a*[-sind(yaw);cosd(yaw);0] , section];
    U1_t = [section - u1_t*[-cosd(yaw)*cosd(az);-sind(yaw)*cosd(az);sind(az)]  , section];
    % U1_t = [U1_a(:,1) - u1_t*[-cosd(yaw)*cosd(az);-sind(yaw)*cosd(az);sind(az)]  , U1_a(:,1)];
    
    plot3(U_inf(1,:),U_inf(2,:),U_inf(3,:),'k','LineWidth',2); hold on;
    plot3(U1_a(1,:),U1_a(2,:),U1_a(3,:),'k--','LineWidth',1); hold on;
    plot3(U1_t(1,:),U1_t(2,:),U1_t(3,:),'k--','LineWidth',1); hold on;
    
end
%% Remove Hold

hold off;