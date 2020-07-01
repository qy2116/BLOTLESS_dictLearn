function num = find_disjoint(SP,eff_pos,I_set,loc)
% find a disjointed column with 'loc' in 'eff_pos'
n = length(eff_pos);
num = [];
for ind = 1:n
    l1 = loc_sparse(SP,I_set,eff_pos(loc));
    l2 = loc_sparse(SP,I_set,eff_pos(ind));
    if isempty(intersect(l1,l2))
        num = ind;
        break
    end
end
    if isempty(num)
        error('More measurements are needed!')
    end
end