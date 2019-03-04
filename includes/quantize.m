function [Xq, indices] = quantize(X, codebook)

    Xq = zeros(size(X));
    indices = zeros(size(X));

    for i = 1:length(X)
        dst = abs(X(i) - codebook);
        [min_val, index] = min(dst);
        Xq(i) = codebook(index);    
        indices(i) = index;
    end

end