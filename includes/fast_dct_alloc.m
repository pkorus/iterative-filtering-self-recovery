function indices = fast_dct_alloc(imgs, block, inblock_mask)
    W = imgs(2);
    H = imgs(1);
    
    S = zeros(H, W);
    
    if size(inblock_mask,1) ~= block
        inb = repmat(inblock_mask, block/size(inblock_mask,1), block/size(inblock_mask,1));
    else
        inb = inblock_mask;
    end

    for bx = 1:W/block
        for by = 1:H/block
            block_id = by + (bx-1)*H/block - 1 + 1;
            S(1+(by-1)*block:(by)*block, 1+(bx-1)*block:(bx)*block) = block_id*64 + inb;
        end
    end
    
    [values, indices] = sort(S(:));
    indices(mod(values,64) == 0) = [];
end
