function [by, bx] = ind2subq(siz, ind)
    by = 1 + mod(ind - 1, siz(1));
    bx = 1 + floor((ind - 1)/siz(1));
end