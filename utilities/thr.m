function x = thr(a)
l = length(a);
for ind = 1:l
    x(ind) = sign(a(ind)) * max(abs(a(ind))-1e-10,0);
end
end