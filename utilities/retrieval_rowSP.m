function rowSP =  retrieval_rowSP(SP)
n = length(SP);
[l,~] = size(SP{1});
rowSP = cell(l,1);
for ind = 1:n
    SP_temp = sum(SP{ind},2);
    for in_ind = 1:l
%         temp = [];
        if SP_temp(in_ind) ==0
            rowSP{in_ind} = [rowSP{in_ind}, ind];
        end
    end
end
end