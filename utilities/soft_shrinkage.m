function x=soft_shrinkage(y,lambda,alpha)
% soft shrinkage operator
  N=length(y);
  x=zeros(N,1);
  for ind=1:N
      x(ind) = sign(y(ind))*max(abs(y(ind))-lambda/alpha,0);
  end
end