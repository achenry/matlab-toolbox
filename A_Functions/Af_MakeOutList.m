function OutList = Af_MakeOutList(IFs)

iList = 1;
%% FAST

OutList{iList} = 'Time'; iList=iList+1;

%% InflowWind
IW_fid     = fopen(([IFs.IW,'.ipt']),'r');

line = '1';
while ~isempty(line)
    line = fgetl(IW_fid);
    if ~ischar(line), break, end
    
    if strfind(line,'OutList')
        
        %start collecting OutList
        while ~isempty(line)
            line = fgetl(IW_fid);
            
            if strfind(line,'END')
                break;
            end
            
            splitLine = strsplit(line);
            OutList{iList} = splitLine{1}; iList=iList+1;
            
            
        end
    end
end

%% ElastoDyn
ED_fid     = fopen(([IFs.ED,'.dat']),'r');

line = '1';
while ~isempty(line)
    line = fgetl(ED_fid);
    if ~ischar(line), break, end
    
    if strfind(line,'OutList')
        
        %start collecting OutList
        while ~isempty(line)
            line = fgetl(ED_fid);
            
            if strfind(line,'END')
                break;
            end
            
            splitLine = strsplit(line);
            OutList{iList} = splitLine{1}; iList=iList+1;
            
            
        end
    end
end

%% AeroDyn
% If using AD14, OutList doesn't appear, so nothing needs to be supressed
% here, no need to include Simulation struct

AD_fid     = fopen(([IFs.AD,'.ipt']),'r');

line = '1';
while ~isempty(line)
    line = fgetl(AD_fid);
    if ~ischar(line), break, end
    
    if strfind(line,'OutList')
        
        %start collecting OutList
        while ~isempty(line)
            line = fgetl(AD_fid);
            
            if strfind(line,'END')
                break;
            end
            
            splitLine = strsplit(line);
            OutList{iList} = splitLine{1}; iList=iList+1;
            
            
        end
    end
end


%% ServoDyn

SD_fid     = fopen(([IFs.SD,'.dat']),'r');

line = '1';
while ~isempty(line)
    line = fgetl(SD_fid);
    if ~ischar(line), break, end
    
    if strfind(line,'OutList')
        
        %start collecting OutList
        while ~isempty(line)
            line = fgetl(SD_fid);
            
            if strfind(line,'END')
                break;
            end
            
            splitLine = strsplit(line);
            OutList{iList} = splitLine{1}; iList=iList+1;
            
            
        end
    end
end

%% Transpose

OutList = OutList';

fclose('all');
