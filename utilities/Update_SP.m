function U = Update_SP(SP,m)
[~,~,l] = size(SP);
U_SP = SP(m+1:end,:,:);
U = cell(1,l);
for ind = 1:l
    d_ind = sum(U_SP(:,:,ind),1)==0;
    temp = U_SP(:,:,ind);
    temp(:,d_ind) = [];
    U{ind} = temp;
end
end