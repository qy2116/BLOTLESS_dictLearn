function D_hat = sort_dic(D)
% Sort dictionary atoms in various descend direction
Var = var(D,1,1);
[~,pos] = sort(Var,'descend');
D_hat = D(:,pos);
end