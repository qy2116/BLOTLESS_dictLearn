function Y = gen_Y_ns(D,X,ind)
len = length(ind);
[m,~] = size(D);
[~,n] = size(X);
Y = zeros(m,n); 
for ii = 1:len
    Y = Y + D(:,ind(ii))*X(ind(ii),:);
end
end