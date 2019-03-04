function bin = binvec2dec(vec)
    ct = 1;
    bin = 0;
    for i=1:length(vec)
        bin = bin + ct*vec(i);
        ct = ct*2;
    end
end