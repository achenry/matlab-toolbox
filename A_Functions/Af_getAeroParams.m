function [Aero,origAerod] = Af_getAeroParams(A_CD,AD15_blade_file)
% Get Aerodynamic parameters from AD15 input blade file


%% Parameters

PLOT = 0;

%% Check file

[~,name,ext] = fileparts(AD15_blade_file);
if isempty(ext)
    AD15_blade_file = [AD15_blade_file,'.dat'];
end

%% Read inFile Parameters
if 0 %DEBUG
    AD15_file = 'SUMR-13_v1_c5_a0333_v113_P139_AeroDyn_blade.dat';
end

fid     = fopen(fullfile(A_CD,'\FAST8_IF\',AD15_blade_file),'r');

line = '1';
while ~isempty(line)
    line = fgetl(fid);
    if ~ischar(line), break, end
    
    
    if strfind(line,'Blade Properties')
        % Numbber of nodes, don't care
        fgetl(fid);
        
        % Property Line
        rawProp = fgetl(fid);
        Properties = strsplit(rawProp);
        
        % Get rid of empty property values
        emptyProperties = cellfun(@isempty,Properties);
        Properties      = Properties(~emptyProperties);
        
        % Unit Line
        rawUnits = fgetl(fid);
        Units   = strsplit(rawUnits);
        
        % Values
        isValue = 1;
        rawValues = '1';
        while true
            rawValues = fgetl(fid);
            if ~ischar(rawValues), break; end
            
            Values{isValue} = strsplit(rawValues);
            isValue = isValue + 1;
            
        end
        break
    end
end

% Sort and shit
for iProp = 1:length(Properties)
    eval([Properties{iProp},'=zeros(1,',num2str(length(Values)),');']);
end

for iNode = 1:length(Values)
    for iProp = 1:length(Properties)
        
        eval([Properties{iProp},'(',num2str(iNode),')=str2double(Values{',num2str(iNode),'}{',num2str(iProp),'});']);
        eval(['Aero.',Properties{iProp},'(',num2str(iNode),')=str2double(Values{',num2str(iNode),'}{',num2str(iProp),'});']);
    end
end



% Compare
if PLOT
    figure(604);
    plot(origBlade.BlFract,[origBlade.BMassDen;Aero.BMassDen],'.-');
    
    figure(605);
    subplot(211);
    plot(origBlade.BlFract,[origBlade.FlpStff;Aero.FlpStff],'.-');
    
    subplot(212);
    plot(origBlade.BlFract,[origBlade.EdgStff;Aero.EdgStff],'.-');
    
    %     figure(606);
    %     subplot(211);
    %     plot(origBlade.BlFract,[origBlade.FlpIner;Blade.FlpIner],'.-');
    %
    %     subplot(212);
    %     plot(origBlade.BlFract,[origBlade.EdgIner;Blade.EdgIner],'.-');
end


%% Close all files
fclose('all');


%% Do Outputs






