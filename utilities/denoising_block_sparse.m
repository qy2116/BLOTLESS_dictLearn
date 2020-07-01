function [psnr_pro, X] = denoising_block_sparse(Y, D, k, ref, sigma)
%%%% Input: Y: noisy image;
%%%%        D: representation dictionary
%%%%        k: block size
%%%% Output: X: denoised image
%%%% Solving a minimization problem: min_{X,z_i,alp_i} = ||Y - X||^2_F
%%%% +sum(||P(X)-D*alp_i||^2_2 + ||alp_i||_1 ) using
%%%% alternating minimization method
%% Initialization
[s1,s2] = size(Y);
[s_d1,s_d] = size(D);
n1 = floor(s1/k);
n2 = floor(s2/k);
alp = zeros(s_d,n1*n2);
X = zeros(n1*k,n2*k);
X_temp = zeros(k*k,n1*n2);
max_iter = 5;
patch_y = patch_img(Y,k);
[patch_y,DC] = removedc(patch_y);
% Y = form_img(patch_y,n1,n2);
% lambda = zeros(s_d1,n1*n2);
for main_iter = 1:max_iter
main_iter

% %% Fix alp, do X update
% Z_w = D*alp;
% % p_y = patch_img(Y,k);
% 
% X_temp = (1*patch_y + 2*sigma/10*Z_w)/(1+2*sigma/10);
% 
% %% Fix X and z, do alp update
% p_x = X_temp;
% % alp = OMPerr(D,p_x,0.009);%30dB
% alp = OMPerr(D,p_x,sigma/1000);%20dB
% % alp = OMPerr(D,p_x,0.5);%10dB
% X = form_img(adddc(p_x,DC),n1,n2);
% psnr_pro(main_iter) = psnr(X,ref);


%% Overlapped Patches
[NN1,NN2] = size(Y);
X_ini = Y;
blocks = im2col(X_ini,[k,k],'sliding');
idx = 1:size(blocks,2);

% go with jumps of 30000
for jj = 1:30000:size(blocks,2)
  
    jumpSize = min(jj+30000-1,size(blocks,2));
    
    %reduceDC
    vecOfMeans = mean(blocks(:,jj:jumpSize));
    blocks(:,jj:jumpSize) = blocks(:,jj:jumpSize) - repmat(vecOfMeans,size(blocks,1),1);
    
    
    Coefs = OMPerr(D,blocks(:,jj:jumpSize),sigma/1000);
   
    %reducedc
    blocks(:,jj:jumpSize)= D*Coefs + ones(size(blocks,1),1) * vecOfMeans;
    
end

count = 1;
Weight = zeros(NN1,NN2);
IMout = zeros(NN1,NN2);
[rows,cols] = ind2sub(size(X_ini)-k+1,idx);
for i  = 1:length(cols)
    col = cols(i); row = rows(i);        
    block =reshape(blocks(:,count),[k,k]);
    IMout(row:row+k-1,col:col+k-1)=IMout(row:row+k-1,col:col+k-1)+block;
    Weight(row:row+k-1,col:col+k-1)=Weight(row:row+k-1,col:col+k-1)+ones(k);
    count = count+1;
end
X_ini = (Y+0.034*sigma*IMout)./(1+0.034*sigma*Weight);

end
X = X_ini;
psnr_pro = 0;
end