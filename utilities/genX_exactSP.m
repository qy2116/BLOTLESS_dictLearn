function X = genX_exactSP(SP)
n = size(SP,2);
[m,~] = size(SP{1});
X = zeros(m,n);
for ind = 1:n
    temp = SP{ind};
    X(:,ind) = SP{ind}*randn(size(temp,2),1);
end
end