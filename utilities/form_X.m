function X = form_X(SP,A)
[m,~,l] = size(SP);
X = zeros(m,l);
for ind = 1:l
    X(:,ind) = SP(:,:,ind)*A(:,ind);
end
end