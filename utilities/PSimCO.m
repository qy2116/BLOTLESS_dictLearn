%Primitive SimCO

function [D] = PSimCO (Y, param)

itN = param.itN;
D = param.initialDictionary;
errglobal = param.errorGoal;
IPara.mu = 0;
IPara.I = param.I;
IPara.dispN = 20; 
IPara.DebugFlag = 0;
IPara.itN = 1;
IPara.gmin = 1e-5; % the minimum value of gradient
IPara.Lmin = 1e-6; %t4-t1 should be larger than Lmin
IPara.t4 = 1e-2; %the initial value of t4
IPara.rNmax = 3; %the number of iterative refinement in Part B in DictLineSearch03.m

for itn = 1:itN

X = OMPerr(D,Y,errglobal);
%X = OMP(D,Y,20);
[D,X] = DictUpdate03 (Y,D,X,IPara);

f=sum(sum((Y-D*X).^2)); 
fprintf('\n ** final cost function value is %3.3e ** \n\n',...
            f);

end

end
