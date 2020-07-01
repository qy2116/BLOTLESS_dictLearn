%Regularized SimCO with mu 0.05 in all 10 iterations

function [D,X] = SimCO(Y, param)

itN = param.iternum;
D = param.initialDictionary;
k = param.Tdata;
% errglobal = param.Edata;
% errglobal = param.errorGoal;
IPara.mu = 0.05;
IPara.I = param.I;
IPara.dispN = 20; 
IPara.DebugFlag = 0;
IPara.itN = 1;
IPara.gmin = 1e-5; % the minimum value of gradient
IPara.Lmin = 1e-6; %t4-t1 should be larger than Lmin
IPara.t4 = 1e-2; %the initial value of t4
IPara.rNmax = 1; %the number of iterative refinement in Part B in DictLineSearch03.m

for itn = 1:itN

    if itn > 5
        IPara.mu = 0.05;
    end
% 
%      if itn < 2
%          IPara.mu = 0.1;
%      end
%      if itn < 4
%          IPara.mu = 0.05;
%      end
%      if itn < 6
%          IPara.mu = 0.03;
%      end    
%      if itn < 8
%          IPara.mu = 0.005;
%      end    

% 
% if itn < 3
%     IPara.mu = 1;
% elseif itn < 5
%         IPara.mu = 0.5;
% elseif itn < 7
%             IPara.mu = 0.1;
%         else 
%             IPara.mu = 0;
% end


% X = OMPerr(D,Y,errglobal);
X = omp(D,Y,[],k);
[D,X] = DictUpdate03 (Y,D,X,IPara);

f=sum(sum((Y-D*X).^2)); 
% fprintf('\n ** final cost function value is %3.3e ** \n\n',f);
end
end

