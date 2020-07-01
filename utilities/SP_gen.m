function SP = SP_gen(m,n)
    % Generate the sparse pattern of DL Test
    % Input: m:dimension of spasce; n: dimension of sparse subspace
    % Output: SP: Sparse Pattern
    index = randperm(m);
    pos = index(1:n);
    SP = zeros(m,n);
    for i = 1:n
        temp = zeros(m,1);
        temp(pos(i)) = 1;
        SP(:,i) = temp;
    end
end