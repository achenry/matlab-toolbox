function [b_op,u_op] = getOpPts(linDir,varargin)
% INPUTS:       linDir  is the directory where linearizations are stored
%                       with prefix Sens 

linFiles = dir(fullfile(linDir,'*.mat'));

for j = 1:length(linFiles)
    
    save_dir = linDir;
    load(fullfile(save_dir,['Sens',num2str(j)]));
    
    u_op(j) = A_OpPt(j,1);
    b_op(j) = str2double(DescCntrlInpt{1}(55:65));
%     S(j)    = AvgDMat(2)*1000;
    
end


% figure(50);
% plot(u_op,rad2deg(b_op)); hold off;
