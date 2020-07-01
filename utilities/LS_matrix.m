function P = LS_matrix(Y,SP,i)
if i==1
    [n,l] = size(Y);
    [m,k] = size(SP{1});
%     [m,k,~] = size(SP);
    P1 = zeros(m*l,m*n);
    P2 = zeros(m*l,k*l);
    %% part1
    for ind = 1:l
        P1(1+m*(ind-1): m*ind,:) = kron(eye(m),Y(:,ind)');
    end

    %% part2
    for ind = 1:l
        temp_vec = zeros(l,1);
        temp_vec(ind) =1;
        P2 = P2 + kron(diag(temp_vec(:,1),0),SP{ind});
    end
    P = [P1 -P2];
else
    [~,l] = size(Y);
    [m,k] = size(SP{1});
    P1 = zeros(m*l,m);
    P2 = zeros(m*l,k*l);
    %% part1
    for ind = 1:l
        P1(m*ind,:) = Y(:,ind)';
    end

    %% part2
    for ind = 1:l
        temp_vec = zeros(l,1);
        temp_vec(ind) =1;
        P2 = P2 + kron(diag(temp_vec(:,1),0),SP{ind});
    end
    P = [P1 -P2];
end