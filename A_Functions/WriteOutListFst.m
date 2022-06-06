% Function: WriteOutListFst FAST-input files
% 
%
% -----------------------------
% Usage:
% -------------
% WriteOutListFst(FSTFile,OutListFst)
% ------------
% Input:
% -------------
% File              String with name of .txt file
% identifier        String of identifier
% value             New value
%
% ------------
% Output:
% ------------
% -
% ------------
% Needs:
% ------------
% -
% ------------
% Modified:
% -------------


% ------------
% ToDo:
% -------------
%
% -----------
% Created: 
% Frank Sandner on 11-Mar-2013
% (c) Universitaet Stuttgart 
% ----------------------------------


function WriteOutListFst(FSTFile, OutList)

TempFSTFile     = [FSTFile(1:end-4), '_temp', FSTFile(end-3:end)]; 
fid             = fopen(FSTFile);
fidTemp         = fopen(TempFSTFile,'w+');
fileending      = FSTFile(end-3:end);

if fid < 0 || fidTemp < 0
    error('WriteOutListFst: Could not open file to be modified or temporary file.')
elseif iscell(OutList) == 0 
    fclose(fid);
    fclose(fidTemp);
    error('WriteOutListFst: OutListFst must have cell format.')
elseif isempty(OutList)
    fclose(fid);
    fclose(fidTemp);
    error('WriteOutListFst: OutListFst is empty, please revise.')
elseif ~strcmp(fileending, '.fst')
    fclose(fid);
    fclose(fidTemp);
    error('WriteOutListFst: File should end on "template.fst", current filename is %s',FSTFile)
else
end

n_list = find(~cellfun(@isempty,OutList));  
n_list = n_list(end);

k_      = [];
linecount = 0;
while isempty(k_) 
    s   = fgetl(fid);
    k_  = strfind(s, 'OutList'); 
    fprintf(fidTemp,'%s\n',s);
    linecount = linecount + 1;
    if linecount > 500
        break;
    end
end

for i_list = 1:n_list
   fprintf(fidTemp, '"%s"\n', OutList{i_list}); 
end
fprintf(fidTemp, 'END of FAST input file (the word "END" must appear in the first 3 columns of this last line).\n');
fprintf(fidTemp, '--------------------------------------------------------------------------------\n');
fprintf(fidTemp, 'Everything below this line has existed in this file before writing the OutList above by WriteOutListFst.m on %s\n', date);
fprintf(fidTemp, '--------------------------------------------------------------------------------\n');

while ~feof(fid)
    s   = fgetl(fid);
    if ~strfind(s, 'END')
        fprintf(fidTemp,'%s\n',s);
    end
end

fclose(fid);
fclose(fidTemp);
recycle = 'on';
delete(FSTFile);
movefile(TempFSTFile,FSTFile);