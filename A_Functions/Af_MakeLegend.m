function hLegend = Af_MakeLegend(introString,sweptParam,outroString,varargin)
% Make legend out of parameter array
% format introString = sweptParam(i)

%% Convert to String
if ~ischar(sweptParam)
    for iParam = 1:length(sweptParam)   %loop over params
        sweptParamString{iParam} = num2str(sweptParam(iParam),4);
    end
end

%% Construct String
for iParam = 1:length(sweptParam)
    legendString{iParam} = [introString,sweptParamString{iParam},outroString];
end

%% Set Legend
if isempty(varargin)
    hLegend = legend(legendString);
else
    hLegend = legend(varargin{1},legendString);
end