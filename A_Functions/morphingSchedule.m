function ConeAngle = morphingSchedule(windSpeed,varargin)

%% Set schedule
U_rated = 11.4;     %m/s
U       = [0,.25*U_rated,.75*U_rated,100*U_rated];
beta    = [2.5,2.5,12.5,12.5];

%% Interpolate
ConeAngle = interp1(U,beta,windSpeed);
