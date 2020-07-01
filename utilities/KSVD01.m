function [D]  = KSVD01(Y,param)

itN = param.itN;
IPara.dispN = 20; 
IPara.DebugFlag = 0;
IPara.gmin=1e-5; 
IPara.itN=1;
IPara.mu = 0;
D = param.initialDictionary;
errglobal = param.errorGoal;

for itn = 1:itN
    
    X = OMPerr(D,Y,errglobal);
    
    [D,X] = K_DictUpdate02 (Y,D,X,IPara);

    f=sum(sum((Y-D*X).^2));
    fprintf('\n ** final cost function value is %3.3e ** \n\n',...
            f);
     
end


end
