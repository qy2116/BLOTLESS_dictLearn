clear all
close all
clc
addpath(genpath('utilities'))
%%
Num_data = 500;                              % Number of samples used                      
load('tdata4.mat');
pos = randperm(15000);
Y = Y_temp(:,pos(1:Num_data));clear Y_temp
Y(:,2:end) = removedc(Y(:,2:end));
Y = column_normalize(Y,2);

%%
flag_ksvd      = 1;
flag_simco     = 1;
flag_blotless  = 1;
flag_mod       = 1;

m = 64;
l = 128;                                     % Dictionary size
k = 10;                                      % Sparsity level
n = Num_data;

%% Initial Dictionary
D_ini = randn(m,l);
D_ini = column_normalize(D_ini,2);

%% Parameters for all algorithms
iter_num = 20;
params.data = Y;
params.Tdata = k;
params.dictsize = l;
params.iternum = iter_num;
params.initdict = D_ini;
params.initialDictionary = D_ini;
params.groudtruthDict = [];
params.odic = D_ini;
params.I= 1:l;
params.verbose = 'no';

%% Dictionary Learning Algorithms
if flag_ksvd
    [Dksvd] = ksvd(params,'');
end

if flag_blotless
    [D_Blot] = Blotless_orig_func(Y, params);
end

if flag_simco    
    [D_simco] = SimCO(Y, params);  
end

if flag_mod   
    D_mod =  MOD_func(Y, params);  
end

%% Save dictionaries
save('Blotless_dict.mat','D_Blot');
save('KSVD_dic.mat','Dksvd');
save('SimCO_dic.mat','D_simco');
save('MOD_dic.mat','D_mod');

%% Show dictionaries
figure;dictshow(D_Blot,'size',[8,8]);title('BLOTLESS')
figure;dictshow(Dksvd,'size',[8,8]);title('K-SVD')
figure;dictshow(D_simco,'size',[8,8]);title('SimCO')
figure;dictshow(D_mod,'size',[8,8]);title('MOD')
