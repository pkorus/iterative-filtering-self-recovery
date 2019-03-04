function vec = dec2binvec(dec, prec)
    if nargin == 1, prec = ceil(log2(dec+1)); end;
    vec = false(1, prec);
    cn = dec;
    i = 1;
    while cn > 0
        bit = mod(cn,2);
        vec(i) = bit;
        cn = floor((cn - bit)/2);
        i = i+1;
    end
    if length(vec) > prec
        vec = vec(1:prec);
        fprintf('WARNING Not enough precision! %d with %d bits\n', dec, prec);
    end
end