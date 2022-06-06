function Af_EditLin(lines, edits, newname, defname, input_mode)
% J. Aho- 3/10/11 jacob.aho@colorado.edu

% This function will take in a vector of .dat file line numbers or a 
% cell array containing strings of Linear parameters 
% (input_mode=1 or 2 respectively) that are to be edited to be edited.
% The defname.dat file is a default Linear file which is copied line by line 
% into newname.dat, replacing lines associated with 'lines' input by the
% associated 'edit' followed by the line of the 'description' reference
% file (not an input, set name below if you want to change this file)

% Description of inputs
% Lines:
% If input_mode==1
%       lines- a vector of line numbers to be edited in X.dat file
% If input_mode==2
%       lines- a cell array w/ each cell as the Linear setting/variable name
% edits- a cell array of strings or numbers- one cell for each element of 'lines' 
% newname- the new Linear file will be called [newname,'.dat']
% defname- the default Linear file which should be called [defname,'.dat']
%    for any line number not an element of (lines), the default file
%    parameters will be copied over.

% Outputs:  Nothing, It makes a new .dat file in the current folder



if (input_mode~=1 && input_mode~=2 )
    error('Error, input_mode must be either 1 or 2')
end

fid=fopen([defname,'.dat']);
if fid==-1
    error(['Error: ', defname, '.dat not found.  Note: you do not need to end string with .dat']);
end

fidW=fopen([newname,'.dat'],'w+');
if fidW==-1
    error(['Error: ', newname, '.dat not found.  Note: you do not need to end string with .dat']);
end

fidD=fopen('FAST_IF\Templates\Linear_desc.dat');
if fidD==-1
    error(['Error: Blank settings AeroDyn input file not found.']);
end

if input_mode==2
    fidI=fopen('FAST_IF\Templates\Linear_inlist.dat');
    if fidI==-1
        error(['Error: Input description file not found.']);
    end
end


if fid>0
    linenum=0;
    editID=1;
    numEdits=length(lines);
    for n=1:numEdits
        if isnumeric(edits{n})
            edits{n}=num2str(edits{n});
        end
    end  
    tline = fgets(fid);
    tlineD = fgets(fidD); %Could optimize this probably
    if input_mode==2
        tlineI=fgets(fidI);
    end
    
    while ischar(tline)
        linenum=linenum+1;
        Fchange=0;
        editID=[];
        if input_mode==1
            editID  = find(lines==linenum, 1, 'last');
            if ~isempty(editID)
                Fchange=1;
                fprintf(fidW,'%s',[edits{editID}, tlineD]);
            end
        elseif input_mode==2
            editID=0;
            for n=1:numEdits
                match=strfind(tlineI,lines{n});
                if match
                    Fchange=1;
                    editID=n;
                    fprintf(fidW,'%s',[edits{editID}, tlineD]);
                end
            end
        end
        if Fchange==0;
            fprintf(fidW,'%s',tline);
        end
        if input_mode==2
        tlineI= fgets(fidI);
        end
        tlineD = fgets(fidD);
        tline = fgets(fid);
    end
end

fclose(fid);
fclose(fidW);
fclose(fidD);
if input_mode==2
    fclose(fidI);
end
disp(['Status: ', newname,'.dat created with ',num2str(numEdits),' changed input(s).']);
