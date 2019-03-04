function bitmask = generate_bitmask(q_map, line_bits)
    bitmask = 2*ones(1, line_bits);
    ci = 1;
    bi = 1;
    while ci <= numel(q_map)

        if q_map(ci) == 0
            blen = 0;
            chk_pattern = [];
        elseif q_map(ci) == 1
            blen = 4;
            chk_pattern = [0 1 1 1];
        elseif q_map(ci) == 2
            blen = 5;
            chk_pattern = [0 0 1 1 1];
        else
            blen = 6;
            chk_pattern = [0 0 1 1 1 1];
        end

        bitmask(bi:bi+blen-1) = chk_pattern;
        bi = bi + blen;
        ci = ci + 1;
    end
end