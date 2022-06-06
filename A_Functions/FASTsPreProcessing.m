%% 1. Initialization 
% clear all;close all;clc;
% run('..\AddingWitlisPaths.m')

%% 2. Loop over processes 
% I only do one process for now, will leave structures in place though
nProcess    = 1;        %FASTProcessingConfig(1);
disp([num2str(nProcess) ' Processes will be preprocessed.'])
for iProcess=1:nProcess
    
    %% 2.1 Get parameters for process 
%     [ProcessName,InitFilesDir,PreProcessingVariation,Disturbance,Parameter]     = FASTProcessingConfig(iProcess);    
    eval(['D',Name_Control(2:end)])
    nInitFiles  = 0;  
    [nVariation, nPermutation, Permutation]     = CreatePermutationMatrix(PreProcessingVariation);
    
    %% 2.2 Loop over permutations 
    for iPermutation    = 1:nPermutation
        
        % 2.2.1 Evaluation PreProcessingVariation 
        for iVariation  = 1:nVariation
            Assignment  = PreProcessingVariation{iVariation,2}(Permutation(iPermutation,iVariation));
            if iscell(Assignment)
                eval([PreProcessingVariation{iVariation,1},'=',   char(Assignment),';'])                
            else
                eval([PreProcessingVariation{iVariation,1},'=',num2str(Assignment),';'])                
            end      
        end                
        
        % 2.2.2 Set up simulation if not already completed
        SimulationNameCell{iPermutation}    = SimulationName;
        if ~exist([InitFilesDir,SimulationName,'_results.mat'],'file');
            
            % Save Init file 
            disp('=======================================================================================');
            fprintf('Writing files for simulation %u out of %u of Process %u ...\n', iPermutation, nPermutation, iProcess);
            save([InitFilesDir,SimulationNameCell{iPermutation},'_init.mat'],'Params');     %Took out Disturbance
            nInitFiles  = nInitFiles+1;
            
%             % Manipulate and copy FAST input files (...not yet)
%             nFASTInputFiles         = length(Parameter.FASTInputFiles);
%             for iFASTInputFiles = 1:nFASTInputFiles                
%                 NewFileName         = Parameter.FASTNewInputFiles{iFASTInputFiles};
%                 TemplateFileName    = Parameter.FASTInputFiles{iFASTInputFiles}; 
%                 disp('-------------------------------------------------------------------------------');
%                 fprintf('FAST template: %s\n', TemplateFileName);
%                 copyfile([Parameter.Simulink.Folder, TemplateFileName] , [InitFilesDir,NewFileName]);
%                 nFASTInputModifications     = size(Parameter.FASTInputModifications{iFASTInputFiles}, 1);
%                 for iFASTInputModifications = 1:nFASTInputModifications
%                     Identifier  = Parameter.FASTInputModifications{iFASTInputFiles}{iFASTInputModifications,2};
%                     NewValue    = Parameter.FASTInputModifications{iFASTInputFiles}{iFASTInputModifications,1};
%                     fprintf('\t\tchanged value of %s into \t%s\n',  Identifier,NewValue);
%                     ManipulateFASTinput([InitFilesDir,NewFileName], Identifier,NewValue);
%                 end
%             end            
%             WriteOutListFst([InitFilesDir,Parameter.FASTNewInputFiles{1}], Parameter.FASTInput.OutList); % Add OutList to fst file
%             
        else
            fprintf('%s already completed, no files written ...\n',SimulationName);
        end
    end

    %% 2.3 save Config for process and clean up 
    save([InitFilesDir,ProcessName,'_config'],'PreProcessingVariation','Permutation','SimulationNameCell')
%     copyfile([cd,'\FASTProcessingConfig.m'],[InitFilesDir,ProcessName,'_FASTProcessingConfig.m']);
    disp([num2str(nInitFiles),' InitFiles written for process ',ProcessName])
    clearvars SimulationNameCell
end


