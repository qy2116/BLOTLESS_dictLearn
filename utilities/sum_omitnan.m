function S = sum_omitnan(R)
[n,~] = size(R);
S = zeros(1,size(R,2));
count = 0;
for ind = 1:n
    if ~isnan(R(ind,:))
        S = S+R(ind,:);
        count = count +1;
    end
end
S = S/count;
end