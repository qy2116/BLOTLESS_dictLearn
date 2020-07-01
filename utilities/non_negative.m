function A_out = non_negative(A)
    mi = min(A);
    A_out = A-mi+0.5;
end