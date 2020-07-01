function [D_ini,err] = Blotless_orig_func(Y,para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Dictionary learning function of BLOTLESS algorithm %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input:                                                     %%%%%%%%%%%
%%%%       Y:    Measured samples                               %%%%%%%%%%%
%%%%       para: A structure contains all the parameters:       %%%%%%%%%%%
%%%%             params.Tdata: column sparisity                 %%%%%%%%%%%
%%%%             params.iternum: total iteration number         %%%%%%%%%%%
%%%%             params.initialDictionary: initial dictionary   %%%%%%%%%%%
%%%%             params.groudtruthDict: ground-truth dictionary %%%%%%%%%%%
%%%%             params.verbose: print the iteration on screen  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Output:                                                    %%%%%%%%%%%
%%%%       D_ini:    Output dictionary                          %%%%%%%%%%%
%%%%       err: Dictionary recovery error at each iteration     %%%%%%%%%%%
%%%%       err is in use only when params.groundtruth is given  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(para.groudtruthDict)
    D = para.groudtruthDict;
    flag_gt = 1;
else
    flag_gt = 0;
end
D_ini = para.initialDictionary;
k = para.Tdata;
[m,l] = size(D_ini);
iter_num = para.iternum;
for out_iter = 1:iter_num
%% Sparse coding 
X_hat = sparse_omp(Y,D_ini,k);
SP = retrieval_SP(X_hat);
num_block = m; % Each block contains m many atoms

for in_iter = 1:(l/num_block)
%% Update S
Ind_s = 1+(in_iter-1)*num_block:num_block+(in_iter-1)*num_block;
Ind_ns = setdiff(1:l,Ind_s);
Y_ns = gen_Y_ns(D_ini,X_hat,Ind_ns);
R = Y - Y_ns;
X_s = X_hat(Ind_s,:);

[s_dim,~] = size(X_s);

for iter_coor = 1:20

LS_mat = [R' X_s'];

if isnan(sum(LS_mat(:)))
    break;
else
   [~,~,v] = svd(LS_mat);
end

V1 = v(1:end-s_dim,end-s_dim+1:end);
V2 = v(end-s_dim+1:end,end-s_dim+1:end);
E = -LS_mat*[V1;V2]*[V1;V2]';
if rank(V2,1e-10)< length(Ind_s)
    break;
end
S = (-V1/(V2))';


%% Update X
lam = 1000;
temp_mat = (LS_mat+E);
X_tp = (S*temp_mat(:,1:m)') - proj_X((S*temp_mat(:,1:m)'),SP,Ind_s,lam,1);

if rank(X_tp,1e-10)< length(Ind_s)
    break;
end
X_ini = X_tp;
X_s = X_ini;

end


%%
D_hat = (temp_mat(:,1:m)')/X_ini;
D_hat = column_normalize(D_hat,2);
D_ini(:,Ind_s) = D_hat;
X_hat(Ind_s,:) = X_ini;
end
if flag_gt
    err(out_iter) = compare_dic(D_ini,D)/l;
end
if  strcmp(para.verbose,'yes') && flag_gt
   fprintf(strcat(sprintf('iter = %i - Dictionary Recovery Error = %2.3f',out_iter, err(out_iter)),'\n'));
end
end
end