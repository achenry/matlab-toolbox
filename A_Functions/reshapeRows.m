%% reshapeRows
% Function: reshapes a Matrix having each row a-times
% 
%% Usage:
%
%  Mout=reshapeRows(M,a)
%
%% Input:
%
% * M- Matrix
% * a- number of multiple rows
%
%% Output:
%
% Mout- output Matrix
%
%% Modified:
%
%
%
%% ToDo:
%
%
%
%% Created: 
% David Schlipf on 2010-02-14
%
% (c) Universitaet Stuttgart
% 


function Mout=reshapeRows(M,a)

sizeM=size(M);
Mout=reshape(repmat(M,1,a)',sizeM(2),sizeM(1)*a)';