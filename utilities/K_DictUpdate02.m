function [D,X,OPara] = K_DictUpdate02(Y,D,X,IPara)

% K_DictUpdate02 is the dictionary update function in K-SVD.
% Given the initial dictionary D, initial sparse coefficient matrix X and
% the traning data matrix Y, this function produces the updated D and X
% by K-SVD algorithm.
%
% References:
% "The K-SVD: An Algorithm for Designing of Overcomplete Dictionaries for 
% Sparse Representation", written by M. Aharon, M. Elad, and A. M. Bruckstein
% and appeared in the IEEE Trans. On Signal Processing, Vol. 54, no. 11, 
% pp. 4311-4322, November 2006
%
%
% See also DictUpdate03
%
% Wei Dai, Tao Xu, Wenwu Wang
% Imperial College London, University of Surrey
% wei.dai1@imperial.ac.uk  t.xu@surrey.ac.uk  w.wang@surrey.ac.uk
% October 2011


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

Omega = X~=0; %record the position of nonzero elements
n = size(Y,2); %the number of training signals
DebugFlag = IPara.DebugFlag;
itN = IPara.itN; 
gmin = IPara.gmin; 
dispN = IPara.dispN; 
OPara.f0 = zeros(1,itN); 
OPara.f1 = zeros(1,itN); 
OPara.gn2 = zeros(1,itN); 
if DebugFlag == 1
    OPara.CondNum = zeros(1,itN); 
end

%% iteration

tic;
[f,X,g] = fg_tilde_eval01(Y,D,Omega,IPara);
for itn = 1:itN
    if itn == 1
        OPara.f0(itn) = f;
    else
        OPara.f0(itn) = OPara.f1(itn-1);
    end
    
    % update dictionary
    jPerm = 1:size(D,2);
    D = findBetterDictionary(Y,D,X,jPerm);
    % save necessary data
    [f,X,g] = fg_tilde_eval01(Y,D,Omega,IPara);
    OPara.f1(itn) = f;
    OPara.gn2(itn) = norm(g,'fro')/norm(Y,'fro')^2;
    if DebugFlag == 1
        Dc=Dcondition(Y,D,Omega);
        OPara.CondNum(itn) = max(Dc(1,:));
    end
    
    % quit condition
    gColn2 = sqrt(sum(g.*g,1));
    gZero = gColn2 < gmin*norm(Y,'fro')^2/n;
    if sum( gZero ) == size(D,2) 
        OPara.f0 = OPara.f0(1:itn);
        OPara.f1 = OPara.f1(1:itn);
        OPara.gn2 = OPara.gn2(1:itn);
        if DebugFlag == 1
            OPara.CondNum = OPara.CondNum(1:itn);
        end        
        return;
    end    
    
    if mod(itn,dispN) == 0
        if DebugFlag == 1
            fprintf('itn=%d: f=%3.3e: gn2=%3.3e: CondNum=%3.3e: time=%f\n',...
                itn,OPara.f1(itn),OPara.gn2(itn),OPara.CondNum(itn),toc);
        else
            fprintf('itn=%d: f=%3.3e: gn2=%3.3e: time=%f\n',...
                itn,OPara.f1(itn),OPara.gn2(itn),toc);
        end
        tic;
    end
end
 
% fprintf('itn=%d: f=%3.3e: gn2=%3.3e:\n',...
%             itn,OPara.f1(itn),OPara.gn2(itn));


function [D,X] = findBetterDictionary(Y,D,X,jPerm)
d = size(D,2);
Omega = X~=0;
for j=jPerm  
    % compute relevant data
    jC = [1:j-1 j+1:d];
    DI = D(:,jC);
    XI = X(jC,:);
    YI = Y - DI*XI;
    
    % update codeword
    OmegaJ = Omega(j,:);
    YIJ = YI(:,OmegaJ);
    [Uj,Sj,Vj] = svds(YIJ,1);
    D(:,j) = Uj;
    X(j,OmegaJ) = Sj*Vj';
end
return;