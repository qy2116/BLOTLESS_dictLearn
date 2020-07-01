function x=least_square(M,beta,z,lam,alpha)
% Solve the least square minimization in the x-subproblem
% INPUT:
%      M: Kernel matrix
%      z: variable z
%      lam: Multiplier
%      beta: linear combination
% OUTPUT:
%      x: The updated estimation of x
% by Qi Yu 22/03/17
[~,num_c] = size(M);
I = eye(num_c);
x=(2*M+alpha*I)\(2*beta-lam'+alpha*z);
end