function Esp = ErrorSP(SP,perc)
Esp = cell(1,size(SP,2));
[m,~] = size(SP{1});
t_num = 0;
B_SP = [];
Err_SP = [];
for ind = 1:size(SP,2)
   temp = SP{ind};
   t_num(ind) = size(temp,2);
   B_SP = [B_SP, temp];
end
sum_num = sum(t_num);
err_num = floor(perc*sum_num);
pos = randperm(sum_num);
pos = pos(1:err_num);

for ind_err = 1:err_num
    sp = eye(m);
    k = randperm(m);
    Err_SP(:,ind_err) = sp(:,k(1));
end

B_SP(:,pos) = Err_SP;

Esp{1} = B_SP(:,1:t_num(1));

for ind = 1:size(SP,2)-1
    st = sum(t_num(1:ind));
    Esp{ind+1} = B_SP(:,1+st:st+t_num(ind));
end

end