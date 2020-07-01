function res = res_loc(loc_square,l)
% take out the loc_sqare from the 1:l set.
ll = 1:l;
res = setdiff(ll,loc_square);
end