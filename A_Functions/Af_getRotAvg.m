function [u_rot,tt,u_rot_filt,MWS] = Af_getRotAvg(inputFileDir,fileName,Parameters)
% Output rotor average wind speed from FF turbulent wind field
global A_CD

if 0 % DEBUG
    
    inputFileDir    = 'D:\sumr13-git-repo\FASTv8\FAST8_IF\WindData\DLC_1_2_14';
    outputFileDir   = 'D:\sumr13-git-repo\FASTv8\SaveData\LidarData';
    
    Parameters.Turbine.String = 'SUMR-13_v1';
    fileName = 'B_10_1.bts';
    
end

load(fullfile('.','TurbineParameters',Parameters.Turbine.String));

%% Read Binary File
if strcmp(inputFileDir(1:2),'F:')
    ADFileName = fullfile(inputFileDir,[fileName,'.bts']);
else
    ADFileName = fullfile('.','FAST8_IF',inputFileDir,[fileName,'.bts']);
end
[velocity, twrVelocity, y, z, zTwr, nz, ny, dz, dy, dt, zHub, z1,MWS]...
    = readTSgrid(ADFileName);


%%
%Streamwise
u_ = velocity(:,1,:,:);
u_bar = reshape(mean(u_,1),[length(y),length(z)]);

%Crosswind
v_ = velocity(:,2,:,:);

if 0
    figure(100);
    imagesc(u_bar');
    colorbar;
    set(gca,'YDir','normal');
end


% single plane
u_plane = squeeze(u_);
v_plane = squeeze(v_);


u_plane_i = squeeze(u_plane(1,:,:));
u_plane_i_bar = mean(u_plane_i(:));


if 0
    figure(101);
    imagesc(u_plane_i')
    title(['Plane Average = ',num2str(u_plane_i_bar,4)])
    colorbar;
    set(gca,'YDir','normal');
end


%% Rotor Plane

zh = Parameters.Tower.Height;

[yy,zz] = meshgrid(y,z);

rotorInd = sqrt(yy.^2 + (zz-zh).^2) < Parameters.Turbine.R * cosd(Parameters.Turbine.ConeAngle);

if 1
    figure(200);
    imagesc(y,z,rotorInd);
end

u_rot_i = mean(u_plane_i(rotorInd));

%% u_rot timeseries
% rotor average is the average over the rotor disk

u_rot = zeros(size(u_plane,1),1);
v_rot = zeros(size(u_plane,1),1);
for tt = 1:size(u_plane,1)
    %streamwise
    u_plane_i   = u_plane(tt,:,:);
    u_rot_i     = mean(u_plane_i(rotorInd));
    u_rot(tt)   = u_rot_i;
    
    %crosswind
    v_plane_i   = v_plane(tt,:,:);
    v_rot_i     = mean(v_plane_i(rotorInd));
    v_rot(tt)   = v_rot_i;
    
    %shear fit
    vert_profile_i = squeeze(mean(u_plane_i,2));
    X = log(z) - log(zh);
    Y = log(vert_profile_i) - log(mean(vert_profile_i));
    alpha_i = X'\Y;
    alpha(tt) = alpha_i;
    
    
    if 0        
        figure(9000);
        plot(X,Y);
        title(num2str(alpha_i));
        keyboard;
    end   
    
end

u_hh = velocity(:,1,ceil(size(velocity,3)/2),ceil(size(velocity,3)/2)); %for comparison
v_hh = velocity(:,2,ceil(size(velocity,3)/2),ceil(size(velocity,3)/2)); %for comparison
tt  = 0:dt:dt*size(velocity,1)-dt;

%% Ideal Filtering

cutoff_freq = .1*2*pi; %rad/s
dt = tt(2) - tt(1);
Nyquist_freq = 1/(2*dt)*2*pi;
n = 500;
wn = cutoff_freq / Nyquist_freq;

% b = fir1(n,wn);
load('b');

filtDelay = (length(b)-1)/2;    %samples

u_rot_filt_ = conv(b,u_rot);
u_rot_filt  = u_rot_filt_(filtDelay+1:end-filtDelay);

v_rot_filt_ = conv(b,v_rot);
v_rot_filt  = v_rot_filt_(filtDelay+1:end-filtDelay);

alpha_filt_ = conv(b,alpha);
alpha_filt  = alpha_filt_(filtDelay+1:end-filtDelay);

if 1
    figure(299);
    subplot(121);
    stem(b);
    
    subplot(122);
    bode(tf(b,1,dt))
    
    figure(300);
    subplot(311);
    plot(tt,u_hh,tt,u_rot);       %looks good, like a smoother value
    hold on;
    plot(tt,u_rot_filt,'LineWidth',2);
    hold off  
    
    
    subplot(312);
    plot(tt,v_hh,tt,v_rot);
    hold on;
    plot(tt,v_rot_filt,'LineWidth',2);
    hold off
    
    subplot(313);
    plot(tt,alpha,tt,alpha_filt);
    hold on;
    plot(tt,alpha_filt,'LineWidth',2);
    hold off
end




%% Save Data
if strcmp(inputFileDir(1:2),'F:')
    save(fullfile(inputFileDir,[fileName,'_lidar']));
else
    save(fullfile('.','FAST8_IF',inputFileDir,[fileName,'_lidar']));
end






