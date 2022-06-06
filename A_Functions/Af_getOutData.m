function [OutData,OutList,OutUnit] = Af_getOutData(fileDir,fileName)
%% Read Output Data
fileName = [fileName,'.out'];
disp(['Reading data from ',fullfile(fileDir,fileName),'...'])
OutData = dlmread(fullfile(fileDir,fileName),'\t',9,0);

%% Read in List and Units

%Open File
fid = fopen(fullfile(fileDir,fileName));

%Make String Formatting
A = size(OutData,2);
a_format = '';
for aa = 1:A
    a_format = [a_format,'%s\t'];
end

%Init Cell Arrays
OutList = cell(1,A);
OutUnit = cell(1,A);
for aa = 1:9
    %Read in a couple lines
    tline = fgets(fid);
    if aa == 7
        temp = textscan(tline,a_format);
        
        %Clean up textscan ouput
        for i = 1:length(temp)
            OutList{i} = temp{i}{1};
        end
    end
    
    if aa == 8
        temp = textscan(tline,a_format);
        
        %Clean up textscan ouput
        for i = 1:length(temp)
            OutUnit{i} = temp{i}{1};
        end
    end
    
end

%Close File
fclose(fid);