function [D] = M_DictUpdate02(Y,X)

% M_DictUpdate02 is the dictionary update function in MOD.
% Given the initial dictionary D, initial sparse coefficient matrix X and
% the traning data matrix Y, this function produces the updated D
% by MOD algorithm.
%
% References:
% "Method of Optimal Directions for Frame Design", written by K. Engan in Proc. ICASSP 1999, 
% pp. 2443-2446, March 1999
%
%
% See also DictUpdate03, K_DictUpdate02
%
% Wei Dai, Tao Xu, Wenwu Wang
% Imperial College London, University of Surrey
% wei.dai1@imperial.ac.uk  t.xu@surrey.ac.uk  w.wang@surrey.ac.uk
% March 2012


%IPara.DebugFlag : the flag for debugging ill-conditioned problem
%   The default value is 0. Set it to 1 for debugging
%IPara.itN : the number of iterations. default value is 400 for synthetic tests 
%IPara.gmin : the minimum value of gradient. default value is 1e-5
%IPara.dispN : the display interval. default value is 20
%OPara.f0 : the initial cost function value
%OPara.f1 : the final cost function value
%OPara.gn2 : the frobenius norm of the gradient
%OPara.CondNum : the condition number of trained dictionary



%% initialization
% 
% Omega = X~=0; %record the position of nonzero elements
% [m,n] = size(Y);
% d = size(D,2);
% itN = IPara.itN; 
% DebugFlag = IPara.DebugFlag;
% gmin = IPara.gmin;
% dispN = IPara.dispN;
% OPara.f0 = zeros(1,itN); 
% OPara.f1 = zeros(1,itN); 
% OPara.gn2 = zeros(1,itN); 
% if DebugFlag == 1
%     OPara.CondNum = zeros(1,itN); 
% end
% 
% %% iteration
% 
% tic;
% 
% %%compute initial f
% X = zeros(d,n);
% for cn = 1:n
%     X(Omega(:,cn),cn) = ...
%         D(:,Omega(:,cn)) \ Y(:,cn);
% end
% Yr = Y - D*X;
% f = sum(sum(Yr.*Yr));
% 
% for itn = 1:itN
%     if itn == 1
%         OPara.f0(itn) = f;
%     else
%         OPara.f0(itn) = OPara.f1(itn-1);
%     end
%     
%     % update dictionary
%     D = findBetterDictionary(Y,X);
%     
%     %% compute updated X f g
%     X = zeros(d,n);
%     for cn = 1:n
%         X(Omega(:,cn),cn) = ...
%         D(:,Omega(:,cn)) \ Y(:,cn);
%     end
%     Yr = Y - D*X;
%     % the cost function 
%     f = sum(sum(Yr.*Yr));
%     g = -2*Yr*X';
%     % additional steps to make sure the orthoganilty
%     DGcorr = sum(D.*g,1);
%     g = g - D.*repmat(DGcorr,m,1);
%     
%     %%output
%     OPara.f1(itn) = f;
%     OPara.gn2(itn) = norm(g,'fro')/norm(Y,'fro')^2;
%     if DebugFlag == 1
%         Dc=Dcondition(Y,D,Omega);
%         OPara.CondNum(itn) = max(Dc(1,:));
%     end
%     
%     % quit condition
%     gColn2 = sqrt(sum(g.*g,1));
%     gZero = gColn2 < gmin*norm(Y,'fro')^2/n;
%     if sum( gZero ) == size(D,2) 
%         OPara.f0 = OPara.f0(1:itn);
%         OPara.f1 = OPara.f1(1:itn);
%         OPara.gn2 = OPara.gn2(1:itn);
%         if DebugFlag == 1
%             OPara.CondNum = OPara.CondNum(1:itn);
%         end        
%         return;
%     end    
%     
%     if mod(itn,dispN) == 0
%         if DebugFlag == 1
%             fprintf('itn=%d: f=%3.3e: gn2=%3.3e: CondNum=%3.3e: time=%f\n',...
%                 itn,OPara.f1(itn),OPara.gn2(itn),OPara.CondNum(itn),toc);
%         else
%             fprintf('itn=%d: f=%3.3e: gn2=%3.3e: time=%f\n',...
%                 itn,OPara.f1(itn),OPara.gn2(itn),toc);
%         end
%         tic;
%     end
% end
%  
% fprintf('itn=%d: f=%3.3e: gn2=%3.3e:\n',...
%             itn,OPara.f1(itn),OPara.gn2(itn));
D = findBetterDictionary(Y,X);
%%subfunction
function D = findBetterDictionary(Y,X)

d = size(X,1); %the number of atoms in dictionary
[m,n] = size(Y); % find dimension of Y

% m = size(Y,1);
% D=Y*X'/(X*X' + 1e-7*speye(size(X,1)));
% SumAtoms = sum(abs(D));
% zerosIdx = find(SumAtoms<eps);
% D(:,zerosIdx) = randn(size(D,1),length(zerosIdx));

%Dictionary = Data*CoefMatrix'*inv(CoefMatrix*CoefMatrix' + 1e-7*speye(size(CoefMatrix,1)));
%sumDictElems = sum(abs(Dictionary));
%zerosIdx = find(sumDictElems<eps);
%Dictionary(:,zerosIdx) = randn(size(Dictionary,1),length(zerosIdx));

% update dictionary
D = (X'\Y')';
% A = X * X';
% B = Y * X';
% D = B / A;

% handle error events: cope with NAN & INF
if sum(sum(isnan(D)))>0 || sum(sum(isinf(D)))>0 
    D = Y(:,ceil(rand(1,d)*n))+0.1*randn(m,d);
    D = D - ones(m,1)*mean(D);
end
em = 0;
temp = sum(isnan(D));
if (sum(temp) > 0)
    em = 1;
else
    temp = sum(isinf(D));
    if (sum(temp) > 0)
        em = 1;
    else
        temp = sum(D.*D);
        if (max(temp) > 10)
            em = 1;
        elseif min(temp) < 0.1
            em = 1;
        end
    end
end


if (em==1)
    D = Y(:,ceil(rand(1,d)*n))+0.1*randn(m,d);
    D = D - ones(m,1)*mean(D);
end

%normalize dictionary
sd = sqrt(sum(D.*D));
D = D.*repmat(1./sd,m,1);
%nsd = 1./sd;
%D = D.*(ones(m,1)*nsd);


return;