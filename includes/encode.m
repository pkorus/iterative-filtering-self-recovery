function imgw = encode(I, line_bits, Qeff, secret_key, settings)

    if exist('codebooks.mat', 'file')
        load codebooks.mat
    else
        throw(MException('Encoder:MissingData', 'File codebooks.mat not found!'));
    end

    if settings.display == 1
        wh = waitbar(0.0, 'Encoding: Setting up protection structures');
    end


    I = im2double(I);
    [h, w] = size(I);

    rand('twister', secret_key);

    % Define representation fidelity for successive coefficients
    HASH_LEN = 8;
    if settings.bit_alloc_metric == 2
        q_map = allocate_line_descriptor_bits_L2(line_bits - HASH_LEN);
    else
        q_map = allocate_line_descriptor_bits(line_bits - HASH_LEN);
    end

    % q_map = [4 3*ones(1,10) 2*ones(1,6) 1*ones(1,6)];
    % bitmask = [zeros(1,2) ones(1,4) repmat([zeros(1,2) ones(1,4)], 1, 10) repmat([zeros(1,2) ones(1,3)], 1, 6) repmat([zeros(1,1) ones(1,3)], 1, 6) zeros(1,HASH_LEN)];

        
    SYN_LENGTH = 6*sum(q_map==4) + 6*sum(q_map==3) + 5*sum(q_map==2) + 4*sum(q_map==1);
    LINES = (numel(I)/32/32);
    % LINES = floor(LINES * 8);


    % Set 1 for bits that need to be hashed - we want to check only MSBs -
    % allow for minor errors
    % bitmask = [ones(1, SYN_LENGTH - HASH_LEN) zeros(1,HASH_LEN)];
    % TODO - NEED TO UPDATE THE BITMASK!
    % bitmask = [zeros(1,2) ones(1,4) repmat([zeros(1,2) ones(1,4)], 1, 10) repmat([zeros(1,2) ones(1,3)], 1, 6) repmat([zeros(1,1) ones(1,3)], 1, 6) zeros(1,HASH_LEN)];
    bitmask = generate_bitmask(q_map, line_bits);

    % Set permutation map - most important bits should be embedded in lower
    % frequences
    permutation = [find(bitmask==2) find(bitmask==1) find(bitmask==0)];
    [~, inv_permutation] = sort(permutation);

    signatures = zeros(LINES, SYN_LENGTH + HASH_LEN, 'uint8');

    fprintf('Line descriptor: %d bits\n', SYN_LENGTH + HASH_LEN);
    fprintf('Line count: %d\n', LINES);

    % -------------- Sample lines and generate descriptors -------------------=

    % f1 = fopen('debug1.dat', 'w');
    lengths = zeros(1,LINES);

    for k = 1:LINES
        
        best_sum = Inf;
    %     best_sum = 0;
        coords_best = zeros(4,1);
        
        % TODO delta in PRNG seed
        for l = 1:1
            coords = rand(4,1);
            x1 = max(1, round(coords(1)*w));
            x2 = max(1, round(coords(2)*w));
            y1 = max(1, round(coords(3)*h));
            y2 = max(1, round(coords(4)*h));

            x = linspace(x1,x2,1000);
            y = linspace(y1,y2,1000);
            index = sub2ind(size(I),round(y),round(x));    
            scan_profile = I(index);
            scan_spectrum = dct(scan_profile);
            
    %         new_sum  = sum(abs(scan_spectrum(2:end)));
            new_sum = sqrt((x1-x2)^2 + (y1-y2)^2);
            if new_sum < best_sum
                coords_best = coords;
                best_sum = new_sum;
            end
        end
        coords = coords_best;
            x1 = max(1, round(coords(1)*w));
            x2 = max(1, round(coords(2)*w));
            y1 = max(1, round(coords(3)*h));
            y2 = max(1, round(coords(4)*h));

            x = linspace(x1,x2,1000);
            y = linspace(y1,y2,1000);
            index = sub2ind(size(I),round(y),round(x));    
            scan_profile = I(index);
            scan_spectrum = dct(scan_profile);
        
            
        % quantize coefficients
        % 1. dc 6 bits uniform
        % 2. ac 6/5/4 bits Lloyd
        scan_spectrumq = zeros(size(scan_spectrum));
        scan_spectrumi = zeros(size(scan_spectrum));
        for i = 1:4
            [scan_spectrumq(q_map == i), scan_spectrumi(q_map == i)] = quantize(scan_spectrum(q_map == i), codebook{i});
        end
        
        % convert to bits
        li = 1;
                    
        % Add coefficients to the signature
        for i = 1:length(q_map)
            if q_map(i) == 4
                len = 6;
            elseif q_map(i) == 3
                len = 6;
            elseif q_map(i) == 2
                len = 5;
            elseif q_map(i) == 1
                len = 4;
            end
            signatures(k, li:(li+len-1)) = dec2binvec(scan_spectrumi(i)-1, len);
            li = li + len;
        end
        
        % Hash
        data = signatures(k, bitmask == 1);    
        crc_data = crc9(double(data));
        len = HASH_LEN;
        signatures(k, li:(li+len-1)) = crc_data;
        signatures(k, :) = signatures(k, permutation);
        
    %     fwrite(f1, sprintf('%s -> %s\n', sprintf('%d', data), sprintf('%d', crc_data)));
        
        % Update image based on the current profile
        line_length = (sqrt((x1-x2)^2 + (y1-y2)^2));
        lengths(k) = line_length;
        
        if settings.display == 1
            waitbar(0.05 + k/LINES*0.75, wh, 'Encoding: preparing line descriptors');
        end

    end

    % fclose(f1);

    % ----------------- Embedding ---------------------------------------------

    if settings.display == 1
        waitbar(0.80, wh, 'Encoding: preparing embedding structures');
    end


    bits_per_block = line_bits/16;
    if bits_per_block ~= round(bits_per_block)
        throw(MException('Encoder:InvalidParams', 'Line precision needs to be a multiple of 16!'));
    end

    fprintf('Using %d bits per block (%.3f bpp)\n', bits_per_block, bits_per_block/64);

    if settings.first_coeff+bits_per_block > 64
        throw(MException('Encoder:InvalidParams', 'Embedding capacity exceeded!'));
    end

    EMB_MASK = zeros(8,8);
    c_mask_ind = zigzag(8);
    EMB_MASK(c_mask_ind(settings.first_coeff:settings.first_coeff+bits_per_block-1)) = 1;
    % EMB_MASK = EMB_MASK';

    % Enumerate components starting from lower frequencies
    c_mask = zeros(8,8);
    c_mask(c_mask_ind) = 1:64;
    EMB_MASK = EMB_MASK.*c_mask;
            
    BLOCK = 32;
            
    T = dctmtx(8);

    Ti = inv(T);
    f = @(x) T*x*T';
    fi = @(x) Ti*x*Ti';
    S = blkproc(I, [8 8], f);

    signatures = signatures';
    bitstream = signatures(:)';
    % signatures = signatures';
            
    fprintf('Bit-stream length: %d bits\n', length(bitstream));

    % select the coefficients for watermarking
    indices = fast_dct_alloc(size(I), BLOCK, EMB_MASK);

    % Verify if there is enough capacity
    fprintf('Allocated capacity: %d\n', length(indices));
    fprintf('Required capacity: %.2f%%\n', 100*length(bitstream)/length(indices));

    if length(indices) ~= length(bitstream)
        throw(MException('Encoder:InvalidParams','Paylaod size and selection channel mismatch!'));
    end

    if settings.display == 1
        waitbar(0.85, wh, 'Encoding: embedding');
    end


    base_quantas = floor(S(indices)/Qeff/2)';
    quantas = [base_quantas-1; base_quantas; base_quantas+1];
    candidates = 2*Qeff*quantas + repmat(Qeff*double(bitstream),3,1);
    distortions = abs(repmat(S(indices)', 3, 1) - candidates);
    [~, mindices] = min(distortions);
    S(indices) = candidates(sub2ind(size(candidates), mindices, 1:size(candidates,2)));     

    imgw = blkproc(S, [8 8], fi);
    imgw = im2uint8(imgw);

    if settings.display == 1
        waitbar(1.0, wh, 'Encoding: protection completed');
        close(wh);
    end

end
