function [a_hat,R,y_hat,error,R_adj,kInd] = paramEstLLS(x,y,varargin)
% with offset
% x is [numObservations x numStates]
% y is [numObservations x 1]
% varargin{1} -> figure number to plot
% varargin{2} -> output n largest outliers

if isempty(varargin)
    PLOT = 0;
else
    PLOT = 1;
end

%varargin
if ~isempty(varargin)
    figNum = varargin{1};
end

%% Error Catching
if size(x,2) > size(x,1)
    x = x';
    disp('Need Tall X...transposing');
end

if size(y,2) > size(y,1);
    y = y';
    disp('Need Tall Y...transposing');
end




%% LLS Param Estimation

X = [x,ones(size(x,1),1)];
a_hat = inv(X'*X)*X'*y;


%% Calculate Residual (R^2)

y_hat = X*a_hat;
error = abs(y - y_hat);

R = 1 - sum(error.^2)/sum(abs((y - mean(y))).^2);

%% Calculate Adjusted Residual 

n = length(y);
d = size(x,2);

nd = (n-1)/(n-d-1);

R_adj = 1 - sum(error.^2)/sum(abs((y - mean(y))).^2) * nd;


%% Plotting
if PLOT
    if exist('figNum','var');
        figure(figNum);
    else
        figure(101);
    end
    
    if size(x,2) == 2
        
        h = plot3(x(:,1),x(:,2),y,'o','MarkerSize',10);
        h.MarkerFaceColor = 'b';
        h.MarkerEdgeColor = 'k';
        grid on;
        
        [xx,yy] = meshgrid(linspace(min(x(:,1)),max(x(:,1)),5),linspace(min(x(:,2)),max(x(:,2)),5));
        
        M_est = a_hat(1)*xx + a_hat(2)*yy + a_hat(3);
        
        hold on;
        hs = surf(xx,yy,M_est,ones(size(x)));
        hs.FaceAlpha = 0.5;
        hs.FaceColor = 'r';
        hs.LineStyle = 'none';
        
        hold off
        
        title(['R = ',num2str(R)]);
    elseif size(x,2) == 1
        plot(x,y,'.');
        
        xx = linspace(min(x),max(x));
        hold on;
        plot(xx,a_hat(1) * xx + a_hat(2));
        hold off;
        title(['R = ',num2str(R)]);
    end
    
end


%% Output Largest Outliers

if length(varargin) == 2
    kOutliers = varargin{2};
    
    [errorSort,indSort] = sort(error);
    
    kInd = indSort(end-kOutliers:end);
end
    
    
    
    
    



