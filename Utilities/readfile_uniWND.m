function [wnd_time, dt, u, u_dir, v, horiz_shear, pwr_law_vert_shear, lin_vert_shear, gust_speed, upflow_angle] = readfile_uniWND(FileName)
%[velocity, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars] = readfile_WND(FileName)
% formerly readBLgrid()
% Input:
% FileName       - string, containing file name to open (.wnd extension is optional)
%
% Output:
%  velocity      - 4-D vector: time, velocity component, iy, iz 
%  y             - 1-D vector: horizontal locations y(iy)
%  z             - 1-D vector: vertical locations z(iz)
%  nz, ny        - scalars: number of points in the vertical/horizontal
%                  direction of the grid
%  dz, dy, dt    - scalars: distance between two points in the vertical [m]/
%                  horizontal [m]/time [s] dimension
% zHub           - hub height [m]
% z1             - vertical location of bottom of grid [m above ground level]
% SumVars        - variables from the summary file (zHub, Clockwise, UBAR, TI_u, TI_v, TI_w)

%%
len    = length(FileName);
ending = FileName(len-3:len);

if strcmpi( ending, '.wnd' )
    FileName = FileName(1:len-4);
end

%-------------------------------------------------------------

 % initialize variables
fileFmt  = 'int16';
ConvFact = 1.0; %results in meters and seconds

str      = {'TIME','WIND SPEED','WIND DIR','VERTICAL SPEED','HORIZ SHEAR','PWR LAW VERT SHEAR', 'LIN VERT SHEAR', 'GUST SPEED', 'UPFLOW ANGLE'};  %MUST be in UPPER case
numVars  = length(str);
data = [];
SummVars = zeros(numVars, 1);

%% -----------------------------------------
%  READ THE HEADER OF THE BINARY FILE 
%  ----------------------------------------- 

wnd_time = [];
u = [];
u_dir = [];
v = [];
horiz_shear = [];
pwr_law_vert_shear =[ ];
lin_vert_shear = [];
gust_speed = []; 
upflow_angle = [];

fid_wnd   = fopen( [ FileName '.wnd' ] );
if ( fid_wnd <= 0 )
   error( 'Wind file could not be opened.' );
end
line = fgetl(fid_wnd); % description header
line = fgetl(fid_wnd); % blank
line = fgetl(fid_wnd); % measurements
line = fgetl(fid_wnd); % measurements
line = fgetl(fid_wnd); % units

k = 1;
while ~feof(fid_wnd)
%     [wnd_time(k), u(k), u_dir(k), v(k), horiz_shear(k), ...
%         pwr_law_vert_shear(k), lin_vert_shear(k), gust_speed(k), ...
%         upflow_angle(k)] = ...
        line = sscanf(fgetl(fid_wnd), '%f', [1, 9]);
        wnd_time(k) = line(1);
        u(k) = line(2);
        u_dir(k) = line(3);
        v(k) = line(4);
        horiz_shear(k) = line(5);
        pwr_law_vert_shear(k) = line(6);
        lin_vert_shear(k) = line(7);
        gust_speed(k) = line(8);
        upflow_angle(k) = line(9);

    k = k + 1;
end
fclose(fid_wnd);
dt = wnd_time(2) - wnd_time(1);
end