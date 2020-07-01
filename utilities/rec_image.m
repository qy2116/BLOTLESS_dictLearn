function im_patch = rec_image(Y)
[m,n] = size(Y);
im_patch = zeros(sqrt(m),sqrt(m),n);
for ind = 1:n
    im_patch(:,:,ind) = reshape(Y(:,ind),sqrt(m),sqrt(m));
end
end