function [OutData] = fastout2OutData(fastout_original)
%FASTOUT2OUTDATA Summary of this function goes here
%   Detailed explanation goes here

OutData.time = fastout_original.Time;

OutData.signals.dimensions = numel(fieldnames(fastout_original))-2;

fields = fieldnames(fastout_original);
fastout_cell = struct2cell(fastout_original);
for i = 3:1:numel(fields)
    
    OutData.signals.values(:,i-2) = fastout_cell{i};
end

end