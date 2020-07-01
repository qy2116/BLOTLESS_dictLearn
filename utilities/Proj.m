function P_y =  Proj(Y,D,SP)
[m,l] = size(Y);
P_y = zeros(m,l);
for ind = 1:l
    pos_temp = SP(1:m,:,ind);
    pos_temp = sum(pos_temp,2);
    pos = pos_temp~=0;
    D_p = D(:,pos);
    if size(D_p,2)==1
        temp = sum((Y(:,ind).'*D_p).*D_p,2);
        P_y(:,ind) = Y(:,ind) - temp;
    elseif isempty(D_p)
        P_y(:,ind) = Y(:,ind);
    else
        D_p = GramSchmidt(D_p);
        temp = sum((Y(:,ind).'*D_p).*D_p,2);
        P_y(:,ind) = Y(:,ind) - temp;
    end
end
end