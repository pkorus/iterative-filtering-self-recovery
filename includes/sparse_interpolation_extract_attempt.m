function [skipped] = sparse_interpolation_extract_attempt(coeffs, Qeffr, LINES, HASH_LEN, secret_key, inv_permutation, bitmask)
    % Extract the lines, verify their integrity, end recover

    res_bitstream = mod(floor(coeffs/Qeffr+0.5),2)';
    signatures_rec = reshape(res_bitstream, [length(res_bitstream)/LINES LINES ])';

    rand('twister', secret_key);

    skipped = 0;

    for k = 1:LINES

        signature = signatures_rec(k, inv_permutation);        
        crc = signature(end-HASH_LEN+1:end); 

        % check crc
        data = signature(bitmask == 1);    
        crc_data = crc9(double(data));

        if ~all(crc == crc_data) 
            skipped = skipped + 1;
            continue;
        end
    end

end