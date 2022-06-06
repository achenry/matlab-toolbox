function [ ] = writeTSgrid( ADFileName, fileFmt, descrStrng, velocity, twrVelocity, nz, ny, dz, dy, dt, zHub, z1, mffws )
% Author: Eric Simley, University of Colorado, Jan. 13, 2014, based on code
% written by Bonnie Jonkman, National Renewable Energy Laboratory.
% This function writes a TurbSim .bts wind file using the wind speeds stored 
% in "velocity". The code performs the inverse of readTSgrid written by 
% Bonnie Jonkman, National Renewable Energy Laboratory. Most of the code is
% adapted from readTSgrid.

% Input:
%  ADFileName    - string: contains file name to create
%  fileFmt       - string: contains format of grid points, must
%                  be a signed int variable: int8, int16, int32, int64
%  descrStrng    - string: optional desciption of wind file
%  velocity      - 4-D array: time, velocity component (1=U, 2=V, 3=W), iy, iz 
%  twrVelocity   - 3-D array: time, velocity component, iz
%  y             - 1-D array: horizontal locations y(iy)
%  z             - 1-D array: vertical locations z(iz)
%  zTwr          - 1-D array: vertical locations of tower points zTwr(iz)
%  nz, ny        - scalars: number of points in the vertical and horizontal
%                  directions of the grid
%  dz, dy, dt    - scalars: distance between two points in the vertical
%                  [m], horizontal [m], and time [s] dimensions
% zHub           - scalar: hub height [m]
% z1             - scalar: vertical location of bottom of grid [m above ground level]
% mffws          - scalar: mean hub-height wind speed

fid  = fopen( ADFileName, 'w' );

tmp = fwrite(fid,7,'int16');

nt = size(velocity,1);

if (isempty(twrVelocity))
    ntwr = 0;
else
    ntwr = size(twrVelocity,3);
end

tmp = fwrite( fid, nz, 'int32');        % the number of grid points vertically, INT(4)
tmp = fwrite( fid, ny, 'int32');        % the number of grid points laterally, INT(4)
tmp = fwrite( fid, ntwr, 'int32');        % the number of tower points, INT(4)
tmp = fwrite( fid, nt, 'int32');        % the number of time steps, INT(4)

tmp = fwrite( fid, dz, 'float32');      % grid spacing in vertical direction, REAL(4), in m
tmp = fwrite( fid, dy, 'float32');      % grid spacing in lateral direction, REAL(4), in m
tmp = fwrite( fid, dt, 'float32');      % grid spacing in delta time, REAL(4), in m/s
tmp = fwrite( fid, mffws, 'float32');      % the mean wind speed at hub height, REAL(4), in m/s
tmp = fwrite( fid, zHub, 'float32');      % height of the hub, REAL(4), in m
tmp = fwrite( fid, z1, 'float32');      % height of the bottom of the grid, REAL(4), in m

%%%calculate Vslope and Voffset variables

if (strcmpi(fileFmt,'int8'))
    intmax = 2^7 - 1;
    intmin = -2^7;
elseif (strcmpi(fileFmt,'int16'))
    intmax = 2^15 - 1;
    intmin = -2^15;    
elseif (strcmpi(fileFmt,'int32'))
    intmax = 2^31 - 1;
    intmin = -2^31;    
elseif (strcmpi(fileFmt,'int64'))
    intmax = 2^63 - 1;
    intmin = -2^63;  
    
else
    intmax = 1;
    intmin = 0;
end

%%%u component
velmax = max(max(max(velocity(:,1,:,:))));
velmin = min(min(min(velocity(:,1,:,:))));

vals = inv([velmax, 1; velmin, 1])*[intmax; intmin];
Vslope(1) = vals(1);
Voffset(1) = vals(2);

%%%v component
velmax = max(max(max(velocity(:,2,:,:))));
velmin = min(min(min(velocity(:,2,:,:))));

vals = inv([velmax, 1; velmin, 1])*[intmax; intmin];
Vslope(2) = vals(1);
Voffset(2) = vals(2);

%%%w component
velmax = max(max(max(velocity(:,3,:,:))));
velmin = min(min(min(velocity(:,3,:,:))));

vals = inv([velmax, 1; velmin, 1])*[intmax; intmin];
Vslope(3) = vals(1);
Voffset(3) = vals(2);

if strcmpi(fileFmt,'float32')
    Voffset = 0.0*Voffset;
    Vslope  = ones(size(Vslope));
end

tmp = fwrite( fid, Vslope(1), 'float32'); % the U-component slope for scaling, REAL(4)
tmp = fwrite( fid, Voffset(1), 'float32'); % the U-component offset for scaling, REAL(4)
tmp = fwrite( fid, Vslope(2), 'float32'); % the V-component slope for scaling, REAL(4)
tmp = fwrite( fid, Voffset(2), 'float32'); % the V-component offset for scaling, REAL(4)
tmp = fwrite( fid, Vslope(3), 'float32'); % the W-component slope for scaling, REAL(4)
tmp = fwrite( fid, Voffset(3), 'float32'); % the W-component offset for scaling, REAL(4)

% Write the description string: e.g. "Generated by TurbSim (vx.xx, dd-mmm-yyyy) on dd-mmm-yyyy at hh:mm:ss."

tmp = fwrite( fid, length(descrStrng), 'int32');     % the number of characters in the description string, max 200, INT(4)
tmp = fwrite( fid,  int8(descrStrng), 'int8' ); % the ASCII integer representation of the character string

% write the velocities

if (isempty(twrVelocity))
    nvTwr = 0;
else
    nvTwr = size(twrVelocity,3);
end

for it = 1:nt


    for iz = 1:nz
        for iy = 1:ny
            for k=1:3
                tmp = fwrite( fid,  Vslope(k)*velocity(it,k,iy,iz) + Voffset(k), fileFmt );
                
            end %k
        end %iy
    end % iz

    %---------------------
    %get the tower points
    %---------------------
    if nvTwr > 0

        for iz=1:nvTwr
            for k=1:3      % scale the data
                tmp = fwrite( fid,  Vslope(k)*twrVelocity(it,k,iz) + Voffset(k), fileFmt );
                
            end   
        end

    end

end %it

fclose(fid);

end

