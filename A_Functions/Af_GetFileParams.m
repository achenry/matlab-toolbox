function [params] = Af_GetFileParams(lines, defname)

%% Open Files
fidDefault  = fopen(defname,'r');


%% Copy Default Lines Over
paramLine = false;
params = cell(1,length(lines));

iLine = 1;
while 1
    tline = fgetl(fidDefault);
    paramLine = false;
    %See if one of lines{} is in default file line
    for iLines = 1:length(lines)
        ntline = length(tline);
        lineInd = strfind(tline,lines{iLines});
        
        if lineInd
            paramLine = true;
            tempLine = strsplit(tline);
            params{iLine} = tempLine{1}; iLine=iLine+1;
        end
        
    end
    
    
    %Break @ end of file
    if ~ischar(tline)
        break;
    end
end


%% Close Files

fclose('all');

