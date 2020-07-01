function X_e = pres_en(X,X_ini)
[~,n] = size(X_ini);
X_e = zeros(size(X_ini));
for ind = 1:n
    fac = norm(X_ini(:,ind),2)/norm(X(:,ind),2);
    X_e(:,ind) = X(:,ind)*fac;
end
end