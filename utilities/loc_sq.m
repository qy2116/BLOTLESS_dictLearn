function loc = loc_sq(SP,I_set,n)
[~,~,l] = size(SP);
pos_ind = zeros(1,l);
for ind = 1:l
    ss = setdiff(1:n,I_set);
    pos_ind(:,ind) =sum(sum(SP(ss,:,ind),1));
end
loc = find(pos_ind==0);
end