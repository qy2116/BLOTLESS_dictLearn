function [A_hat,col_norm] = column_normalize(A,k)
    % Function Description
    % Column normalize of a real matrix A
    %
    % Input: A: Matrix
    % Output: A_hat: column_normalized A using k-norm
    % Created by Qi Yu 21/02/2018
    % Last modified 06/03/2018
    % log: modified to k-norm column_normalize, use vecnorm(.) to replace norm(.)
    [m1,m2] = size(A);
    A_hat = zeros(m1,m2);
    col_norm = zeros(m2,1);
    for ind = 1:m2
        col_norm(ind) = norm(A(:,ind),k);
        A_hat(:,ind) = A(:,ind)./col_norm(ind);
    end
end