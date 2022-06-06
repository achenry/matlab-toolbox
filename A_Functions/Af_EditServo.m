function Af_EditServo(lines, edits, newname, defname, templateDir)

%% Open Files
fidDefault  = fopen([defname,'.dat'],'r');
fidNew      = fopen([newname,'.dat'],'w+');
fidDesc     = fopen(fullfile(templateDir,'ServoDyn_desc.dat'));

%% Make All Edits into string
for iEdits = 1:length(edits)
    if ~ischar(edits)
        edits{iEdits} = num2str(edits{iEdits});
    end
end

%% If Edits are empty, copy file over and close
if isempty(edits)
    copyfile([defname,'.dat'],[newname,'.dat']);
else
    %% Copy Default Lines Over
    editedLine = false;
    while 1
        tline = fgetl(fidDefault);
        editedLine = false;
        %See if one of lines{} is in default file line
        for iLines = 1:length(lines)
            if ~isempty(strfind(tline,lines{iLines}))
                if strfind(tline,lines{iLines}) < 35
                    editedLine = true;
                    %find description
                    fidDesc     = fopen(fullfile(templateDir,'ServoDyn_desc.dat'));
                    while 1
                        tlineDesc   = fgetl(fidDesc);
                        if ~isempty(strfind(tlineDesc,lines{iLines}))
                            %write new value and description
                            fprintf(fidNew,'%s%s\n',edits{iLines},tlineDesc);
                            break;
                        end
                        if ~ischar(tlineDesc)
                            break;
                        end
                    end
                    fclose(fidDesc);
                end
            end
            
        end
        
        %Break @ end of file
        if ~ischar(tline)
            break;
        end
        
        if ~editedLine
            fprintf(fidNew,'%s\n',tline);
        end
        
        
        
    end
    
end

%% Close Files

fclose('all');
disp(['Status: ', newname,'.dat created with ',num2str(length(edits)),' changed input(s).']);

