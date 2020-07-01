function pos = find_eff(SP,I_set,p)
% find the active set contains the m+p row in sparse pattern SP
    pos = [];
%     [mm,~,~] = size(SP);
    temp = SP(p,:,:);
    s_temp = sum(temp);
    pp = setdiff(I_set,p);
    temp2 = SP(pp,:,:);
    s_temp2 = sum(sum(temp2,1));
    pos = intersect(find(s_temp~=0),find(s_temp2==0));
    if isempty(pos)
        error('More measurements are needed.')
    end
end