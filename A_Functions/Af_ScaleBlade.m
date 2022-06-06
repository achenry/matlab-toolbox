function [outFile,Blade,origBlade] = Af_ScaleBlade(eta,inFile,Blade)
% Scale blade using aeroelatic scaling
% eta is the length scaling factor
% inFile is the original ElastoDyn blade input
% Blade is the blade struct from Parameters.Blade
% Blade.k_M is the mass scaling (0 for linear mass scaling, 1 for cubic
% mass scaling)
% Blade.k_Fs is the flapwise stiffness scaling (0 for linear stiffness scaling, 1
% for 5th power)
% Blade.k_Es is the edgewise stiffness scaling (0 for linear stiffness scaling, 1
% for 5th power)
% Inertia scale factors will be the same as the mass factors for now, i.e.
% k_FI = k_M, k_EI = k_M
%% Parameters

PLOT = 0;

if ~isfield(Blade,'k_M')
    Blade.k_M = 0;
end

if ~isfield(Blade,'k_Fs')
    Blade.k_Fs = 0;
end

if ~isfield(Blade,'k_Es')
    Blade.k_Es = 0;
end
%% Read inFile Parameters
if 0 %DEBUG
    inFile = 'SUMR-13i_v6_s4_Blade.dat';
end

fid     = fopen(fullfile('.','\FAST8_IF\',inFile),'r');

outFile = [inFile(1:7),'_M',num2str(Blade.k_M),'_Fs',num2str(Blade.k_Fs),'_Es',num2str(Blade.k_Es),'.dat'];
fid2    = fopen(fullfile('.','.\FAST8_IF\',outFile),'w');

line = '1';
while ~isempty(line)
    line = fgetl(fid);
    if ~ischar(line), break, end
    
    if strfind(line,'FAST INDIVIDUAL BLADE FILE')
        fprintf(fid2,'%s\n',line);
        fgetl(fid); %advance one line and print new header
        fprintf(fid2,'%s\n',['Base file: ',inFile,', eta = ',num2str(eta),', k_M = ',num2str(Blade.k_M),...
            ', Fs = ',num2str(Blade.k_Fs),', Es = ',num2str(Blade.k_Es)]);
        
    else
        fprintf(fid2,'%s\n',line);
    end
    
    if strfind(line,'DISTRIBUTED BLADE PROPERTIES')
        % Property Line
        rawProp = fgetl(fid);
        Properties = strsplit(rawProp);
        fprintf(fid2,'%s\n',rawProp);
        
        % Unit Line
        rawUnits = fgetl(fid);
        Units   = strsplit(rawUnits);
        fprintf(fid2,'%s\n',rawUnits);
        
        % Values
        isValue = 1;
        while isValue
            rawValues = fgetl(fid);
            if strfind(rawValues,'BLADE MODE SHAPES')
                isValue = 0;
            else
                Values{isValue} = strsplit(rawValues);
                isValue = isValue + 1;
            end
            
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
        eval(['origBlade.',Properties{iProp},'(',num2str(iNode),')=str2double(Values{',num2str(iNode),'}{',num2str(iProp),'});']);
    end
end

%% Scaling Parameters

Del     = eta - 1;
k_M     = Blade.k_M;
k_Fs    = Blade.k_Fs;
k_Es    = Blade.k_Es;
k_FI    = k_M;      %Hardcoded like this
k_EI    = k_M;      %ditto

MassScaleExp    = 2;
% MassScaleFactor = k_M*(1+Del)^MassScaleExp + (1-k_M)
MassScaleFactor = eta^(MassScaleExp*k_M);

% k_Fs            = 1;
% k_Es            = 1;
StiffScaleExp   = 4;
% FlapStiffFactor = k_Fs*(1+Del)^StiffScaleExp + (1-k_Fs);
% EdgeStiffFactor = k_Es*(1+Del)^StiffScaleExp + (1-k_Es);
FlapStiffFactor = eta^(StiffScaleExp*k_Fs);
EdgeStiffFactor = eta^(StiffScaleExp*k_Es);


k_FI            = k_M;
k_EI            = k_M;
InertScaleExp   = 4;
% FlapInerFactor = k_FI*(1+Del)^StiffScaleExp + (1-k_FI);
% EdgeInerFactor = k_EI*(1+Del)^StiffScaleExp + (1-k_EI);
FlapInerFactor = eta^(InertScaleExp*k_M);
EdgeInerFactor = eta^(InertScaleExp*k_M);


Blade.BlFract    = origBlade.BlFract;
Blade.BMassDen   = origBlade.BMassDen * MassScaleFactor;
Blade.FlpStff    = origBlade.FlpStff * FlapStiffFactor;
Blade.EdgStff    = origBlade.EdgStff * EdgeStiffFactor;

if isfield(origBlade,'AeroCent')
    Blade.AeroCent   = origBlade.AeroCent;
    Blade.StrcTwst   = origBlade.StrcTwst;
    Blade.GJStff     = origBlade.GJStff * eta^4;
    Blade.EAStff     = origBlade.EAStff * eta^2;
    Blade.Alpha      = origBlade.Alpha;
    Blade.FlpIner    = origBlade.FlpIner * FlapInerFactor;
    Blade.EdgIner    = origBlade.EdgIner * EdgeInerFactor;
    Blade.PrecrvRef  = origBlade.PrecrvRef;
    Blade.PreswpRef  = origBlade.PreswpRef;
    Blade.FlpcgOf    = origBlade.FlpcgOf;
    Blade.EdgcgOf    = origBlade.EdgcgOf;
    Blade.FlpEAOf    = origBlade.FlpEAOf;
    Blade.EdgEAOf    = origBlade.EdgEAOf;
end

% Compare
if PLOT
    figure(604);
    semilogy(origBlade.BlFract,[origBlade.BMassDen;Blade.BMassDen],'.-');
    
    figure(605);
    subplot(211);
    semilogy(origBlade.BlFract,[origBlade.FlpStff;Blade.FlpStff],'.-');
    
    subplot(212);
    semilogy(origBlade.BlFract,[origBlade.EdgStff;Blade.EdgStff],'.-');
    
    %     figure(606);
    %     subplot(211);
    %     plot(origBlade.BlFract,[origBlade.FlpIner;Blade.FlpIner],'.-');
    %
    %     subplot(212);
    %     plot(origBlade.BlFract,[origBlade.EdgIner;Blade.EdgIner],'.-');
end


%%  Write new blade parameters
if Blade.k_M ~= 0 || Blade.k_Fs ~= 0 || Blade.k_Es ~= 0
    
    if isfield(origBlade,'AeroCent')  %long format or short format?
        % Write properties
        for iNode = 1:length(Values)
            fprintf(fid2,'%1.5f\t%1.3f\t%1.3f\t%1.3f\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t\n',...
                [Blade.BlFract(iNode),Blade.AeroCent(iNode),Blade.StrcTwst(iNode),Blade.BMassDen(iNode),...
                Blade.FlpStff(iNode),Blade.EdgStff(iNode),Blade.GJStff(iNode),Blade.EAStff(iNode),...
                Blade.Alpha(iNode),Blade.FlpIner(iNode),Blade.EdgIner(iNode),Blade.PrecrvRef(iNode),...
                Blade.PreswpRef(iNode),Blade.FlpcgOf(iNode),Blade.EdgcgOf(iNode),Blade.FlpEAOf(iNode),...
                Blade.EdgEAOf(iNode)]);
        end
        
    else
        for iNode = 1:length(Values)
            fprintf(fid2,'%1.5f\t%1.3f\t%1.3f\t%1.3f\t%1.5e\t%1.5e\n',...
                [Blade.BlFract(iNode),origBlade.PitchAxis(iNode),origBlade.StrcTwst(iNode),Blade.BMassDen(iNode),...
                Blade.FlpStff(iNode),Blade.EdgStff(iNode)]);
        end
        
        
        
    end
    
    fprintf(fid2,'%s\n',rawValues);  %(Blade mode shape header)
    
    while ischar(line)
        line = fgetl(fid);
        fprintf(fid2,'%s\n',line);
    end
    
else
    disp(['Blade Not Scaled, keeping same name']);
    outFile = inFile;
end
%% Close all files
fclose('all');


%% Do Outputs






