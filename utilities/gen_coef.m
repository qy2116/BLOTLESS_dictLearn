function X = gen_coef(l,n,k)
X = zeros(l,n);
for ind = 1:n
    pos = randperm(l);
    c_pos = k;
    temp = zeros(1,l);
%     temp(pos(1:c_pos(1))) = spli_gaussian(1,c_pos(1),0.5);
    temp(pos(1:c_pos(1))) = randn(1,c_pos(1));
    X(:,ind) = temp;
end
end
