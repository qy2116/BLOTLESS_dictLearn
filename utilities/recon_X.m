function X = recon_X(SP)
n_SP = size(SP,2);
[m,~] = size(SP{1});
X = zeros(m,n_SP);
for ind = 1:n_SP
    X(:,ind) = SP{ind}*ones(size(SP{ind},2),1);
end
end