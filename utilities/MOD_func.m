function [D]  = MOD_func(Y,param)

itN = param.iternum;
k = param.Tdata;
D = param.initialDictionary;
% errglobal = param.Edata;%param.errorGoal;

for itn = 1:itN
    
%     X = OMPerr(D,Y,errglobal);
%     X = OMPerr(D,Y,errglobal);%OMPerr_qi(D, Y, errglobal);
    X = omp(D,Y,[],k);
    D = M_DictUpdate02 (Y,X);
    f=sum(sum((Y-D*X).^2));
%     fprintf('\n ** final cost function value is %3.3e ** \n\n',f);
     
end

end
