function B_I = image_block(I,s,num)
[m,n] = size(I);
B_I = zeros(s*s,num);
for ind = 1:num
    rand_pos1 = randperm(m-s+1);
    rand_pos2 = randperm(n-s+1);
    pos1 = rand_pos1(1);
    pos2 = rand_pos2(1);
    B(:,:,ind) = imcrop(I,[pos1 pos2 s-1 s-1]);
    B_I(:,ind) = reshape(B(:,:,ind),s*s,1);
end
end