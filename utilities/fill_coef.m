function X = fill_coef(SP)
n = length(SP);
[l,k] = size(SP{1});
X = zeros(l,n);
for ind = 1:n
    temp = randn(k,1);
    X(:,ind) = SP{ind}*temp;
end
end