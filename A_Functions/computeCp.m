function Cp = computeCp(lambda,beta,a_hat)

Cp = zeros(size(lambda));
for i = 1:size(lambda,2)
    x = [lambda(:,i),    lambda(:,i).^2,  beta(:,i),   beta(:,i).^2,...
        lambda(:,i) .* beta(:,i),    lambda(:,i) .* beta(:,i).^2,...
        lambda(:,i).^2 .* beta(:,i),  lambda(:,i).^2 .* beta(:,i).^2,   ones(length(lambda),1)];
    Cp(:,i) =  x * a_hat;
end

% force positive values of Cp?
thresh = 0.025;
Cp(Cp<thresh) = thresh;

maxThresh = 0.5;
Cp(Cp>maxThresh) = maxThresh;
% for i = 1:nl
%     y = lambda.^i;
% end
%
% for i = 1:nb
%     y = y + P.^i;
% end
% k = nl+nb;
% for i = 1:nl
%     for j = 1:nb
%
%         X(:,k+1) = L.^i .* P.^j;
%         k = k+1;
%     end
% end