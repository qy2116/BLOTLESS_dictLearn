function I = isallone(a)
a = abs(a);
l = length(a);
for ind =1:l
   a(ind) = a(ind)-1;
   a(ind) = max(abs(a(ind))-1e-8,0);
end
    I = isequal(a,zeros(l,1));
end