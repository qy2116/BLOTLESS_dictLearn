function P = get_nz_pos(X)
[l,n] = size(X);
P = zeros(l,n);
for ind = 1:l
    P(ind,(X(ind,:)~=0)) = 1; 
end
end