function SP_block = find_SP(SP,Ind_s)
n = length(SP);
SP_block = cell(1,n);
for ii = 1:n
   tt = SP{ii};
   temp = tt(Ind_s,:);
   SP_block{ii} = temp(:,(sum(temp,1)~=0)); 
end



% pos = [];
% for ii = 1:n
%     sp = SP{ii};
%     sum_sp = sum(sp,2); 
%     if sum(sum_sp(Ind_ns)) == 0
%         pos = [pos ii];
%     end
% end
% 
% if isempty(pos)
%     error('More measurements are needed!')
% else
%     ll = length(pos);
%     SP_block = zeros(l,k,ll);
%     for ind = 1:ll
%         SP_block(:,:,ind) = SP{pos(ind)};
%     end
% end
end