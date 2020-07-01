function R = adddc(I,DC)
[~,n] = size(I);
for ind = 1:n
    R(:,ind) = I(:,ind) + DC(ind);    
end
end