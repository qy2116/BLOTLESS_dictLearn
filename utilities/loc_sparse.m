function loc = loc_sparse(SP, I_set, pos)
%locate the sparse positions in SP for a specific pos
temp = SP(I_set,:,pos);
tt = sum(temp,2);
loc = find(tt~=0);
% loc = setdiff(loc,nsq_set);
end