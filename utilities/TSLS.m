function D = TSLS(Y,param)
D_ini = param.initialDictionary;
[m,n] = size(Y);
Iternum = param.itN;
k = 8;
l = param.K;
errglobal = param.errorGoal;
for out_iter = 1:Iternum
%     out_iter
%% Sparse coding
    
X_hat = sparse_omp(Y,D_ini,k);
s_xhat = sum(X_hat,2);
for reg_ind = 1:length(s_xhat)
    if s_xhat(reg_ind) == 0
        X_hat(reg_ind,:) = 0.1*gen_coef(1,n,1);
    end
end
% X_hat = OMPerr(D_ini,Y,errglobal);
X_hat = full(X_hat);
% SP = retrieval_nSP(X_hat);
SP = retrieval_SP(X_hat);
num_block = m;
for in_iter = 1:(l/num_block)
% out_iter
%% Sparse coding
    
X_hat = sparse_omp(Y,D_ini,k);
%% Add some non-zeros if needed
s_xhat = sum(X_hat,2);
for reg_ind = 1:length(s_xhat)
    if s_xhat(reg_ind) == 0
        X_hat(reg_ind,:) = 0.1*gen_coef(1,n,1);
    end
end


SP = retrieval_SP(X_hat);
% num_block = l/2;
num_block = m;
for in_iter = 1:ceil(l/num_block)
%% Update S
if num_block+(in_iter-1)*num_block < l
    Ind_s = 1+(in_iter-1)*num_block:num_block+(in_iter-1)*num_block;
else
    Ind_s = 1+(in_iter-1)*num_block:l;
end
Ind_ns = setdiff(1:l,Ind_s);
Y_ns = gen_Y_ns(D_ini,X_hat,Ind_ns);
R = Y - Y_ns;
X_s = X_hat(Ind_s,:);
[s_dim,~] = size(X_s);

for iter_coor = 1:10
LS_mat = [R' X_s'];
if isnan(sum(LS_mat(:)))
    continue
else
   [~,sig,v] = svd(LS_mat);
end

V1 = v(1:end-s_dim,end-s_dim+1:end);
V2 = v(end-s_dim+1:end,end-s_dim+1:end);
E = -LS_mat*[V1;V2]*[V1;V2]';
S = (-V1/(V2))';

%% Update X
lam = 1000;
temp_mat = (LS_mat+E);
X_ini = (S*temp_mat(:,1:m)') - proj_X((S*temp_mat(:,1:m)'),SP,Ind_s,lam,1);
X_s = X_ini;
end
%%
D_hat = R/X_ini;
D_hat = column_normalize(D_hat,2);
D_ini(:,Ind_s) = D_hat;
X_hat(Ind_s,:) = X_ini;

end

err(out_iter) = norm(Y-D_ini*X_hat);
end
D = D_ini;
end