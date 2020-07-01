function num = loc_disjoint(SP_res,I_set,eff_pos,sp_loc)
temp = SP_res(:,:,eff_pos);
[~,~,l] = size(temp);
temp = sum(temp,2);
for ind = 1 : l
    if temp(I_set(sp_loc),:,ind) == 0
        num = ind;
        break;
    else
        num = ind;
    end
end
end