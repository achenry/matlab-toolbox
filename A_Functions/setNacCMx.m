function NacCmx = setNacCMx(Parameters,Simulation,ParamString)
global A_CD

defaultNacCmX = Parameters.Turbine.NacCMx;
%% Load Value
% if no ParamString, keep default from sheet
% if no NacelleCenter in SaveData/, keep default from sheet

% if ~isempty(ParamString)
    try
        load(fullfile(A_CD,'SaveData',Parameters.Turbine.String,...
            [Simulation.Name_Control,ParamString],'NacCM_Center','NacelleCenter'));
        useDefault = 0;
        
    catch
        disp(['WARNING setNacCMx(): No NacelleCenter file,  ',num2str(defaultNacCmX),'m for ',...
        Parameters.Turbine.String,' NacCmx']); 
        useDefault = 1;
    end   
% else
%     disp(['WARNING setNacCMx(): ParamString is empty, using ',num2str(defaultNacCmX),'m for ',...
%         Parameters.Turbine.String,' NacCmx']); 
%     useDefault = 1;
% end


%% Set Param

if useDefault
    NacCmx = defaultNacCmX;
else
    NacCmx = xcm_est;
    disp(['Setting NacCmx to ', num2str(NacCmx), ' m, from ',Parameters.Turbine.String,...
            '/',Simulation.Name_Control,ParamString]);
end



