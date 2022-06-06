function Af_EditAero(lines, edits, newname, defname, templateDir, input_mode)
% J. Aho- 3/10/11 jacob.aho@colorado.edu

% This function will take in a vector of .ipt file line numbers or a 
% cell array containing strings of AeroDyn parameters 
% (input_mode=1 or 2 respectively) that are to be edited to be edited.
% The defname.ipt file is a default AeroDyn file which is copied line by line 
% into newname.ipt, replacing lines associated with 'lines' input by the
% associated 'edit' followed by the line of the 'description' reference
% file (not an input, set name below if you want to change this file)

% Description of inputs
% Lines:
% If input_mode==1
%       lines- a vector of line numbers to be edited in X.ipt file
% If input_mode==2
%       lines- a cell array w/ each cell as the AeroDyn setting/variable name
% edits- a cell array of strings or numbers- one cell for each element of 'lines' 
% newname- the new AeroDyn file will be called [newname,'.ipt']
% defname- the default AeroDyn file which should be called [defname,'.ipt']
%    for any line number not an element of (lines), the default file
%    parameters will be copied over.

% Outputs:  Nothing, It makes a new .ipt file in the current folder



if (input_mode~=1 && input_mode~=2 )
    error('Error, input_mode must be either 1 or 2')
end

fid=fopen([defname,'.ipt']);
if fid==-1
    error(['Error: ', defname, '.ipt not found.  Note: you do not need to end string with .ipt']);
end

fidW=fopen([newname,'.ipt'],'w+');
if fidW==-1
    error(['Error: ', newname, '.ipt not found.  Note: you do not need to end string with .ipt']);
end

fidD=fopen(fullfile(templateDir,'AeroDyn_desc.ipt'));
if fidD==-1
    error(['Error: Blank settings AeroDyn input file not found.']);
end

if input_mode==2
    fidI=fopen(fullfile(templateDir,'AeroDyn_inlist.ipt'));
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
disp(['Status: ', newname,'.ipt created with ',num2str(numEdits),' changed input(s).']);
