function [count,pos] = compare_dic(D,D_hat)
count = 0;
[~,l] = size(D);
for ii = 1:l
    for ind1 = 1:l
        dist(ind1) = 1 - abs(D_hat(:,ii)'*D(:,ind1));
    end
    count = count+min(dist);
    [val(ii),pos(ii)] = min(dist);
%     if min(dist) < 0.01
%         count = count+1;
%     end
end
end