function dCpdL = computePartialCp(lambda,beta,a_hat)

dCpdL = zeros(size(lambda));
for i = 1:size(lambda,2)
    x = [ones(length(lambda),1),    2*lambda(:,i),  zeros(length(lambda),1),   zeros(length(lambda),1),...
        beta(:,i),  beta(:,i).^2,...
        2*lambda(:,i).* beta(:,i),  2*lambda(:,i) .* beta(:,i).^2,   zeros(length(lambda),1)];
    dCpdL(:,i) = x * a_hat;
end

