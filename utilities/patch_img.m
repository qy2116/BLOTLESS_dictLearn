function P = patch_img(I,k)
[m,n] = size(I);
m1 = floor(m/k);
m2 = floor(n/k);
num = m1*m2;
P = zeros(k*k,num);
for ind1 = 1:m1
    for ind2 = 1:m2
        temp = I(1+(ind1-1)*k:k+(ind1-1)*k,1+(ind2-1)*k:k+(ind2-1)*k);
        P(:,m2*(ind1-1)+ind2) = reshape(temp,k*k,1);
    end
end
end