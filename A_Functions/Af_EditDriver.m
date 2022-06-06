function Af_EditDriver(lines, edits, newname, defname, templateDir)

%% Open Files
fidDefault  = fopen([defname,'.dvr'],'r');
fidNew      = fopen([newname,'.dvr'],'w+');
fidDesc     = fopen(fullfile(templateDir,'AeroDyn15_desc.dat'));

%% Make All Edits into string
for iEdits = 1:length(edits)
    if ~ischar(edits)
        edits{iEdits} = num2str(edits{iEdits});
    end
end

%% If Edits are empty, copy file over and close
if isempty(edits)
    copyfile([defname,'.dvr'],[newname,'.dvr']);
else
    %% Copy Default Lines Over
    editedLine = false;
    while 1
        tline = fgetl(fidDefault);
        editedLine = false;
        %See if one of lines{} is in default file line
        for iLines = 1:length(lines)
            if ~isempty(strfind(tline,lines{iLines}))
                editedLine = true;
                %find description
                fidDesc     = fopen(fullfile(templateDir,'AeroDyn15_desc.dvr'));
                while 1
                    tlineDesc   = fgetl(fidDesc);
                    if ~isempty(strfind(tlineDesc,lines{iLines}))
                        %write new value and description
                        fprintf(fidNew,'%s%s\n',edits{iLines},tlineDesc);
                    end
                    if ~ischar(tlineDesc)
                        break;
                    end
                end
                fclose(fidDesc);
            end
            
        end
        
        if ~editedLine
            fprintf(fidNew,'%s\n',tline);
        end
        
        
        %Break @ end of file
        if ~ischar(tline)
            break;
        end
    end
    
end

%% Close Files

fclose('all');
disp(['Status: ', newname,'.dvr created with ',num2str(length(edits)),' changed input(s).']);

