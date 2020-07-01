function I = form_img(P,n1,n2)
% Form image from square patches
[m,~] = size(P);
k = sqrt(m);
for ind1 = 1:n1
    for ind2 = 1:n2
        I(1+(ind1-1)*k:k+(ind1-1)*k,1+(ind2-1)*k:k+(ind2-1)*k) = reshape(P(:,n2*(ind1-1)+ind2),k,k);% + DC(n2*(ind1-1)+ind2)*ones(k,k);
    end
end
end