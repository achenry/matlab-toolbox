function ConeAngle = morphingSchedule4(windSpeed,varargin)

%% Set schedule
U_rated = 11.4;     %m/s
U       = [0,   8,      9.5,    100*U_rated];
beta    = [2.5, 2.5,    12.5,   12.5];

%% Interpolate
ConeAngle = interp1(U,beta,windSpeed);
