function Af_EditAero15(lines, edits, newname, defname, templateDir)

%% Open Files
fidDefault  = fopen([defname,'.dat'],'r');
if strcmp(templateDir,'.\AeroDyn_IF\Templates')
    fidNew      = fopen([newname,'.ipt'],'w+');
else
    fidNew      = fopen([newname,'.dat'],'w+');
end
fidDesc     = fopen(fullfile(templateDir,'AeroDyn15_desc.dat'));

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
            ntline = length(tline);
            if ntline > 50
                lineInd = strfind(tline(1:50),lines{iLines});
            else
                lineInd = strfind(tline,lines{iLines});
            end
            
            if lineInd
                editedLine = true;
                %find description
                fidDesc     = fopen(fullfile(templateDir,'AeroDyn15_desc.dat'));
                while 1
                    tlineDesc   = fgetl(fidDesc);
                    ntlineDesc  = length(tlineDesc);
                    if ntlineDesc > 50
                        lineDescInd = strfind(tlineDesc(1:50),lines{iLines});
                    else
                        lineDescInd = strfind(tlineDesc,lines{iLines});
                    end
                    
                    if lineDescInd
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
disp(['Status: ', newname,'.dat created with ',num2str(length(edits)),' changed input(s).']);

