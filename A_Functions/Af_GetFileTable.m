%% Af_getFileTable
% Loads input file templates from Google Sheet "Turbine Model Summary"
%
%
%% Load Sheet
A1_Initialize

GoogleSheetKey = '1vhV8ztQ8W4ToF5pAEhV3PDZC9jCa55UM';
gid_TradeStudy = {'1762862680'};
% {base, SUMR-13, Main,  trade studies}

ModelString = {};
ModelNumber = 0;
for iSheet = 1:length(gid_TradeStudy)
    
    result = GetGoogleSpreadsheet(GoogleSheetKey,gid_TradeStudy{iSheet});
    SheetModels = result(1,2:end);
    
    ModelString = [ModelString , SheetModels];
    FileString = result(:,1);
    
    for iModel = 1:length(SheetModels)
        FAST_Template{iModel+ModelNumber}   = result{strcmp('Fast Input',FileString),strcmp(SheetModels(iModel),result(1,:))};
        ED_Template{iModel+ModelNumber}     = result{strcmp('ElastoDyn Input',FileString),strcmp(SheetModels(iModel),result(1,:))};
        ED_BladeFile{iModel+ModelNumber}    = result{strcmp('Blade Input (ED)',FileString),strcmp(SheetModels(iModel),result(1,:))};
        ED_TowerFile{iModel+ModelNumber}    = result{strcmp('Tower Input',FileString),strcmp(SheetModels(iModel),result(1,:))};
        AD14_Template{iModel+ModelNumber}   = result{strcmp('AeroDyn Input (v14)',FileString),strcmp(SheetModels(iModel),result(1,:))};
        AD15_Template{iModel+ModelNumber}   = result{strcmp('AeroDyn Input (v15)',FileString),strcmp(SheetModels(iModel),result(1,:))};
        AD15_BladeFile{iModel+ModelNumber}  = result{strcmp('Blade Geometry',FileString),strcmp(SheetModels(iModel),result(1,:))};
        SD_Template{iModel+ModelNumber}     = result{strcmp('ServoDyn Input',FileString),strcmp(SheetModels(iModel),result(1,:))};
        IW_Template{iModel+ModelNumber}     = result{strcmp('InflowWind Input',FileString),strcmp(SheetModels(iModel),result(1,:))};
        BD_Template{iModel+ModelNumber}     = result{strcmp('BeamDyn Input',FileString),strcmp(SheetModels(iModel),result(1,:))};
        MD_Template{iModel+ModelNumber}     = result{strcmp('MoorDyn Input',FileString),strcmp(SheetModels(iModel),result(1,:))};
        HD_Template{iModel+ModelNumber}     = result{strcmp('HydroDyn Input',FileString),strcmp(SheetModels(iModel),result(1,:))};

%         if sum(strcmp('HydroDyn Input',FileString))
%             
%         else
%             HD_Template{iModel+ModelNumber}     = '';
%         end
         
        if sum(strcmp('SubDyn Input',FileString))
            SbD_Template{iModel+ModelNumber}     = result{strcmp('SubDyn Input',FileString),strcmp(SheetModels(iModel),result(1,:))};
        else
            SbD_Template{iModel+ModelNumber}     = '';
        end
    end
    
    ModelNumber = ModelNumber + length(SheetModels);
    
    
    
end

%% Save

% if exist('HD_Template','var')
if exist('SbD_Template','var')
    save(fullfile(A_CD,'FAST8_IF','TemplateTable'),'ModelString','FAST_Template',...
        'ED_Template','ED_BladeFile','ED_TowerFile','AD14_Template','AD15_Template','AD15_BladeFile',...
        'SD_Template','IW_Template','BD_Template','HD_Template','MD_Template','SbD_Template');
else
    save(fullfile(A_CD,'FAST8_IF','TemplateTable'),'ModelString','FAST_Template',...
        'ED_Template','ED_BladeFile','ED_TowerFile','AD14_Template','AD15_Template','AD15_BladeFile',...
        'SD_Template','IW_Template','BD_Template','HD_Template','MD_Template');
end
% else
%     save(fullfile(A_CD,'FAST8_IF','TemplateTable'),'ModelString','FAST_Template',...
%         'ED_Template','ED_BladeFile','ED_TowerFile','AD14_Template','AD15_Template','AD15_BladeFile',...
%         'SD_Template','IW_Template','BD_Template');
% end
