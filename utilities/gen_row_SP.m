function row_sp = gen_row_SP(X,i)
% Generate the row sparse pattern of the i th row
[~,n] = size(X);

pos = find(X(i,:)~=0);
l = length(pos);
row_sp = zeros(n,l);
for ind = 1:l
    row_sp(pos(ind),ind) = 1;
end
end 