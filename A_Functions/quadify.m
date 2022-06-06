function [X,X1] = quadify(x1,x2)
% input, set of x values
% output, set of quadratic fitting values (X), with ones (X1), length(X) =
% length(x1) * length(x2)

%% make sure x1, x2 are long and the same length
err = 0;
if size(x1,2) > 1
    x1 = x1';
end

if size(x2,2) > 1
    x2 = x2';
%     disp('Warning in quadify');
end

%% Meshgrid it
[A,B] = meshgrid(x1,x2);

%% Quadify using method 1

if ~err
    X = [A(:),A(:).^2,B(:),B(:).^2,...
        A(:) .* B(:)];
    
    X1 = [X,ones(size(A(:)))];
else
    disp(['quadify() error, check inputs']);
    X = nan;
    X1 = nan;
end