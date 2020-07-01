function [Y_hat,DC] = removedc(Y)
[~,n] = size(Y);
for ind = 1:n
    DC(ind) = mean(Y(:,ind));
    Y_hat(:,ind) = Y(:,ind) - mean(Y(:,ind));    
end
end