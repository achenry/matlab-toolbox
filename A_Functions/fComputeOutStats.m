function R=fComputeOutStats(folder, p, Cases, MatFileName);
% For each .outb file in `folder`, computes statistics of the different channels, including `AdditionalChannels`, and potentially save to `MatFileName`.
%
% 
% INPUTS:
%   folder: location where .outb files will be looked for
%   p: Structure of metrics inputs as returned by fMetricsVars.m
%
% 
% OUTPUTS:
%    The statistics computed and returned in the structure `R` are:
%      - Mean, Max, Std of each channel
%      - Trvl: the "travel" of a signal per second, the time average of the cumulative sum of the absolute increments of the signal.
%      - Spec: the spectra of each channels
% 
% --- Optional arguments
if ~exist('MatFileName','var'); MatFileName = []; end

% --- Reading all output files and computing stats of channels
if ~exist(folder,'file'); error('Folder not found: %s',folder); end;
if folder(end)~='/' && folder(end)~='\'
   folder=[folder '/'];
end
OutFiles = arrayfun(@(x) x.name, dir([folder '*.outb']),'UniformOutput',false);
nFiles  = length(OutFiles);
nCases  = length(Cases.Names);
if nFiles<=0; error('No *.outb files found in folder %s',folder); end;
if nFiles~=nCases; warning('Inconsistent number of files: %d cases are specified in the case file but %d .outb files are found in %s .',nCases,nFiles,folder); end;

% for simplicity, we define the same frequency vector for all files
R.Freq = 0:0.01:2;

R.OutFiles = cell(nCases,1);
for iFile = 1:nCases
    filename=dir([folder Cases.Names{iFile} '*.outb']);
    if isempty(filename);
        error('No file with pattern `%s` were found in folder %s. \n Has this simulation run properly?',[Cases.Names{iFile} '*.outb'],folder);
    elseif length(filename)>1
        error('Several file with pattern `%s` were found in folder %s. Please remove the unwanted files.\n',[Cases.Names{iFile} '*.outb'],folder);
    end
    filename=filename.name;
    R.OutFiles{iFile}=filename;
    fprintf('%s\n', filename);
    [Channels, R.ChanName ] = fReadFASTAddChannels([folder filename], p.AdditionalChannels);
    % --- Computing stats for all channels
    if iFile==1 % allocation
        nChan = length(R.ChanName);
        R.Max  = zeros(nFiles, nChan);
        R.Std  = zeros(nFiles, nChan);
        R.Mean = zeros(nFiles, nChan);
        R.Trvl = zeros(nFiles, nChan);
        R.Spec = cell (nFiles, nChan);
    end
    if nChan~=length(R.ChanName); error('Inconsitent number of channels between files! %d/%d, file: %s in folder %s',nChan,length(R.ChanName),filename,folder); end

    % Time info and selection of data after tStart
    Time = Channels(:,1);
    T  = Time(end)-Time(1);
    dt = Time(2)-Time(1);
    if abs(T-Cases.tSim(iFile))>dt; error('Inconsistent simulation length: %.1f specificid in the case file, but %.1f found in file %s of folder %s',Cases.tSim(iFile),T,filename,folder); end

    [~,iStart] = min(abs(Time-p.tStart));
    ISelect = iStart:size(Channels,1);
    % Stats for all channels 
    for iChan = 1:nChan
        % To save some time now, we only do that for the required channels
        id=find(ismember(p.Vars(:,3),R.ChanName{iChan}));
        if isempty(id)
            continue % We skip this channel
        end
        sig = Channels(ISelect,iChan);
        n=length(sig);
        % --- FFT
        [S,f] =  fpwelch(sig,p.FFT_WinLen,[],[],1/dt,'detrend',true) ; % smoothen more with lower window sizes
        S0 = interp1(f,S,R.Freq);
        % --- Spectral amplitudes at frequencies of interest
        %fprintf('df=%.4f\n',f(2)-f(1))
        % nRefFreq=size(p.FreqVars,1);
        % SpecA=NaN(nRefFreq,1); 
        % for iFreq=1:nRefFreq
        %     SpecA(iFreq,1)=interp1(R.Freq,S0,p.FreqVars{iFreq,2});
        % end
        % figure(iChan)
        % semilogy(f,S,'Color',[0.5,0.5,0.5]); hold all
        % %semilogy(cell2mat(p.FreqVars(:,2)),SpecA,'ko')
        % %semilogy(cell2mat(p.FreqVars(end-4:end,2)),SpecA(end-4:end),'ko','MarkerFaceColor','k')
        % hold all
        % xlim([0,2])
        % title(R.ChanName{iChan})
        % -- Capping power
        if isequal(lower(R.ChanName{iChan}),'genpwr')
            sig(sig>7000)=7000;
        end
        % --- Storing stats and spectrum
        R.Spec{iFile, iChan} = S0;
        R.Mean(iFile, iChan) = mean(sig)            ;
        R.Max (iFile, iChan) = max (sig)            ;
        R.Std (iFile, iChan) = std (sig)            ;
        R.Trvl(iFile, iChan) = sum(abs(diff(sig)))/T;
    end
end

% --- Saving to a matfile
if ~isempty(MatFileName)
    save(MatFileName, '-struct','R');
end



