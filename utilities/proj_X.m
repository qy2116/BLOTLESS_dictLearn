function X_p = proj_X(X_ini, SP, Ind_s, lam, meth)
%%%% when method = 0, soft threshold with known sparse pattern
%%%% when method ~=0, hard threshold with known sparse pattern
if meth == 0
[m,n] = size(X_ini);
X_p = zeros(m,n);
% lam = 1000;
for ind  = 1:n
    sp = SP{ind}'*SP{ind};
    K = eye(m)+lam*sp;
%     K = SP{ind};
    X_p(:,ind) = inv(K)*X_ini(:,ind);
end
else
[m,n] = size(X_ini);
X_p = zeros(m,n);
% lam = 1000;
for ind  = 1:n
    sp = sum(SP{ind},2);
    sp = sp(Ind_s,:);
    sp = ones(m,1)-sp;
    X_p(:,ind) = X_ini(:,ind).*sp;
end 
end
end