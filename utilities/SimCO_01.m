%Regularized SimCO with mu 0.05 in all 10 iterations

function [D,count] = SimCO_01 (Y, param)

itN = param.itN;
D = param.initialDictionary;
sp = param.errorGoal;
D_0 = param.odic;
IPara.mu = 0.05;
IPara.I = param.I;
IPara.dispN = 20; 
IPara.DebugFlag = 0;
IPara.itN = 1;
IPara.gmin = 1e-5; % the minimum value of gradient
IPara.Lmin = 1e-6; %t4-t1 should be larger than Lmin
IPara.t4 = 1e-2; %the initial value of t4
IPara.rNmax = 3; %the number of iterative refinement in Part B in DictLineSearch03.m

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
X = sparse_omp(Y,D,sp);
%X = OMP(D,Y,20);
[D,X] = DictUpdate03 (Y,D,X,IPara);


count(itn) = compare_dic(D,D_0);
% if (count(itn)/length(IPara.I))<0.02
%     break;
% end

f=sum(sum((Y-D*X).^2)); 
% fprintf('\n ** final cost function value is %3.3e ** \n\n',...
%             f);
end
end

