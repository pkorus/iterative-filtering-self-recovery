function [Y, error_map, Qeffr] = decode(Iwm, line_bits, secret_key, settings)

    if exist('codebooks.mat', 'file')
        load codebooks.mat
    else
        throw(MException('Decoder:MissingData', 'File codebooks.mat not found!'));
    end

    if settings.display == 1
        wh = waitbar(0.0, 'Decoding: Setting up protection structures');
    end

    % ---------- Restore params from the encoder ------------------------------

    rand('twister', secret_key);

    [h, w] = size(Iwm);

    % Define representation fidelity for successive coefficients
    HASH_LEN = 8;
    if settings.bit_alloc_metric == 2
        q_map = allocate_line_descriptor_bits_L2(line_bits - HASH_LEN);
    else
        q_map = allocate_line_descriptor_bits(line_bits - HASH_LEN);
    end

    % q_map = [4 3*ones(1,10) 2*ones(1,6) 1*ones(1,6)];
        
    SYN_LENGTH = 6*sum(q_map==4) + 6*sum(q_map==3) + 5*sum(q_map==2) + 4*sum(q_map==1);
    LINES = (numel(Iwm)/32/32);
    bitmask = generate_bitmask(q_map, line_bits);

    permutation = [find(bitmask==2) find(bitmask==1) find(bitmask==0)];
    [~, inv_permutation] = sort(permutation);

    bits_per_block = line_bits/16;
    if bits_per_block ~= round(bits_per_block)
        throw(MException('Decoder:InvalidParams', 'Line precision needs to be a multiple of 16!'));
    end

    if settings.first_coeff+bits_per_block > 64
        throw(MException('Decoder:InvalidParams', 'Embedding capacity exceeded!'));
    end

    BLOCK = 32;
    EMB_MASK = zeros(8,8);
    c_mask_ind = zigzag(8);
    EMB_MASK(c_mask_ind(settings.first_coeff:settings.first_coeff+bits_per_block-1)) = 1;
    % EMB_MASK = EMB_MASK';

    c_mask = zeros(8,8);
    c_mask(c_mask_ind) = 1:64;
    EMB_MASK = EMB_MASK.*c_mask;

    % -------------------------------------------------------------------------

    T = dctmtx(8);

    Ti = inv(T);
    f = @(x) T*x*T';
    fi = @(x) Ti*x*Ti';

    Sd = blkproc(im2double(Iwm), [8 8], f);

    rec_indices = fast_dct_alloc(size(Iwm), BLOCK, EMB_MASK);

    % prepare for watermark retrieval

    % Estimate Quantization ---------------------------------------------------

    % Qeffr = recover_quantization_step_alt(Sd(rec_indices));
    % fprintf('Estimated Q step %.4f\n', Qeffr);

    % Try to fine-tune --------------------------------------------------------

    if settings.display == 1
        waitbar(0.01, wh, 'Decoding: detecting quantization step');
    end


    % if Qeffr == 0
        qes = 0.035:0.001:0.085;
    % else
    %     qes = round(Qeffr*0.5*1000)/1000:0.001:round(Qeffr*2*1000)/1000;
    % end

    errors = zeros(size(qes));
    for qe = qes
        [skipped] = sparse_interpolation_extract_attempt(Sd(rec_indices), qe, LINES, HASH_LEN, secret_key, inv_permutation, bitmask);
        errors(qe == qes) = skipped;
    end

    [~, qind] = min(errors);
    % Qest = Qeffr;
    Qeffr = qes(qind);
    fprintf('Brute-selected Q step %.4f\n', Qeffr);

    % -------------------------------------------------------------------------

    if settings.display == 1
        waitbar(0.03, wh, 'Decoding: extracting watermark');
    end

    res_bitstream = mod(floor(Sd(rec_indices)/Qeffr+0.5),2)';
    signatures_rec = reshape(res_bitstream, [length(res_bitstream)/LINES LINES ])';

    % Extract the lines, verify their integrity -------------------------------

    S = zeros(size(Iwm));
    M = zeros(size(Iwm));

    % rand('twister', secret_key);

    skipped = 0;
    gof_skipped = 0;

    % f1 = fopen('debug2.dat', 'w');

    error_map = zeros(1, LINES);
    
    matlab_version = regexpi(version('-release'), '[0-9]{4}', 'match');
    matlab_version = str2num(matlab_version{1});
    if exist('svmclassify','file') == 2 && matlab_version < 2018
        svm_fallback = false;
    else
        fprintf('WARNING Matlab`s SVM routines not found. Fallback to linear classification.\n');
        svm_fallback = true;
    end


    for k = 1:LINES
        
        signature = signatures_rec(k, inv_permutation);

        coords = rand(4,1);
        x1 = max(1, round(coords(1)*w));
        x2 = max(1, round(coords(2)*w));
        y1 = max(1, round(coords(3)*h));
        y2 = max(1, round(coords(4)*h));
            
        x = linspace(x1,x2,1000);
        y = linspace(y1,y2,1000);
        index = sub2ind(size(Iwm),round(y),round(x));
        
        scan_profile_dct = zeros(1,1000);
        scan_profile_ind = zeros(1,1000);

        % Add coefficients from the signature
        pi = 1;
        li = 1; 
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
            idx = 1+binvec2dec(signature(li:(li+len-1)));
            cbk = codebook{q_map(i)};
            if idx > 0
                scan_profile_dct(pi) = cbk(idx);
            end
            scan_profile_ind(pi) = idx;
            li = li + len;
            pi = pi + 1;
        end
                
        len = HASH_LEN;
        crc = signature(li:(li+len-1)); 
        
        % check crc
        data = signature(bitmask == 1);    
        crc_data = crc9(double(data));

    %     fwrite(f1, sprintf('%s -> %s (recovered %s)\n', sprintf('%d', data), sprintf('%d', crc_data), sprintf('%d', crc)));
        
        if ~all(crc == crc_data) 
            skipped = skipped + 1;
            error_map(k) = 1;
            continue;
        end
        
        mind = find(scan_profile_dct ~= 0, 1, 'last');
        polyf = polyfit([2:mind]', abs(scan_profile_dct(2:mind))', 1);
        
        var = std(abs(scan_profile_dct(2:mind)))^2;
        maxVal = max(abs(scan_profile_dct(2:mind))-polyval(polyf,1:mind-1));
        
        features = [polyf var maxVal];
        
        if settings.display == 1
            waitbar(0.05 + k/LINES*0.5, wh, 'Decoding: validating line descriptors');
        end

        matlab_version = regexpi(version('-release'), '[0-9]{4}', 'match');
        matlab_version = str2num(matlab_version{1});
        if ~svm_fallback
            if svmclassify(svmstr, features) == 0
                gof_skipped = gof_skipped + 1;
                error_map(k) = 2;
                continue
            end
        else
            % Fallback to linear classification if SVM not supported            
            a = [-7/0.4 (6-7/0.4*0.3)];
            distg = (-a(1)*polyf(2)+polyf(1)-a(2))/(sqrt(a(1)^2+a(2)^2));
            if distg < 0
                gof_skipped = gof_skipped + 1;
                continue
            end
        end
        
        scan_profile_compq = idct(scan_profile_dct);

        % Update image based on the current profile
        line_length = (sqrt((x1-x2)^2 + (y1-y2)^2));
        weight = (1/line_length).^10;
        weights = M(index);

        if settings.fusion_type > 1

            [index_u, val_indices, ~] = unique(index);

            for z = 1:length(index_u)
                iu = index_u(z);
                vu = scan_profile_compq(val_indices(z));
                pixel_candidates{iu} = [pixel_candidates{iu} vu];
                pixel_weights{iu} = [pixel_weights{iu} weight];
            end

            M(index) = 1;

        else

            if settings.fusion_type == 0
                S(index) = S(index).*(weights./(weights+weight)) + scan_profile_compq.*(weight./(weights+weight));    
                M(index) = M(index) + weight;   % put weight here            
            else                
                replace_at = M(index) < weight;
                S(index(replace_at)) = scan_profile_compq(replace_at);
                M(index(replace_at)) = weight;
            end

        end
    end

    % fclose(f1);

    fprintf('using %d lines (%d crc-skipped + %d fit-skipped)\n', LINES-skipped-gof_skipped, skipped, gof_skipped);

    if skipped+gof_skipped > 0.995*LINES
        if settings.display == 1
            waitbar(1.0, wh, 'Decoding: decoding failed');
            close(wh);
        end
        
    %     throw(MException('Decoder:AuthError', 'More than 99.5%% of the lines skipped, stopping!'));
        Y = uint8(127*ones(size(Iwm)));
        error_map = reshape(error_map, size(Y)/32);
        return;
    end

    if settings.fusion_type > 1
        
        if settings.display == 1
            waitbar(0.55, wh, 'Decoding: fusing candidate values');
        end
        
        for iu = 1:numel(S)
            candidates = pixel_candidates{iu};
            if isempty(candidates) 
                continue
            end
            if length(candidates) == 1
                S(iu) = candidates;
                continue
            end
            weights = pixel_weights{iu};
            weights = weights/sum(weights);
            if settings.fusion_type == 2
                S(iu) = mean(candidates);
            elseif settings.fusion_type == 3
                S(iu) = sum(candidates.*weights);
            end
        end
    end
        
    M = double(M ~= 0);

    S(S < 0) = 0;

    % ------------------ Begin Content Reconstruction -------------------------

    if settings.disable_borders ~= 1

        border_dmax = 3;

        [im_h, im_w] = size(S);

        % Sample each border using 32 points - put value from closes known point
        x_sam = round(linspace(1, im_w, 64));
        y_sam = round(linspace(1, im_h, 64));

        % Find closest known point for UPPER boundary -----------------------------

        if settings.display == 1
            waitbar(0.6, wh, 'Decoding: Predicting borders');
        end

        ci = 1;
        Mtemp = M;
        Mtemp(ci+1:end,:) = 0;

        while sum(Mtemp(:)) < im_h/4
            ci = ci + 1;
            Mtemp = M;
            Mtemp(ci+1:end,:) = 0;
        end

        [known_y, known_x] = find(Mtemp ~= 0);
        y2 = 1;

        for x2 = x_sam
            best_distance = Inf;
            pixel_value = NaN;
            for i = 1:length(known_x)
                y1 = known_y(i);
                x1 = known_x(i);
                distance = sqrt((x1-x2)^2 + (y1-y2)^2);
                if distance < best_distance
                    best_distance = distance;
                    pixel_value = S(y1, x1);
                end
            end
            S(y2,x2) = pixel_value;
        end

        samples = S(y2, x_sam);
        for d = 0:border_dmax
            S(y2+d,:) = spline(x_sam, samples, 1:im_w);
            M(y2+d,:) = 1;
        end

        % Find closest known point for BOTTOM boundary ----------------------------

        if settings.display == 1
            waitbar(0.625, wh, 'Decoding: Predicting borders');
        end

        
        ci = 1;
        Mtemp = M;
        Mtemp(1:end-ci,:) = 0;

        while sum(Mtemp(:)) < im_h/4
            ci = ci + 1;
            Mtemp = M;
            Mtemp(1:end-ci,:) = 0;
        end

        [known_y, known_x] = find(Mtemp ~= 0);
        y2 = im_h;

        for x2 = x_sam
            % find closest known point for bottom boundary
            best_distance = Inf;
            pixel_value = NaN;
            for i = 1:length(known_x)
                y1 = known_y(i);
                x1 = known_x(i);
                distance = sqrt((x1-x2)^2 + (y1-y2)^2);
                if distance < best_distance
                    best_distance = distance;
                    pixel_value = S(y1, x1);
                end
            end
            S(y2,x2) = pixel_value;
        end

        samples = S(y2, x_sam);
        for d = 0:border_dmax
            S(y2-d,:) = spline(x_sam, samples, 1:im_w);
            M(y2-d,:) = 1;
        end

        % Find closest known point for RIGHT boundary -----------------------------

        if settings.display == 1
            waitbar(0.65, wh, 'Decoding: Predicting borders');
        end

        ci = 1;
        Mtemp = M;
        Mtemp(:,1:end-ci) = 0;

        while sum(Mtemp(:)) < im_w/4
            ci = ci + 1;
            Mtemp = M;
            Mtemp(:,1:end-ci) = 0;
        end

        [known_y, known_x] = find(Mtemp ~= 0);
        x2 = im_w;

        for y2 = y_sam
            % find closest known point for bottom boundary
            best_distance = Inf;
            pixel_value = NaN;
            for i = 1:length(known_x)
                y1 = known_y(i);
                x1 = known_x(i);
                distance = sqrt((x1-x2)^2 + (y1-y2)^2);
                if distance < best_distance
                    best_distance = distance;
                    pixel_value = S(y1, x1);
                end
            end
            S(y2,x2) = pixel_value;
        end

        samples = S(y_sam, x2);
        for d = 0:border_dmax
            S(:,x2-d) = spline(y_sam, samples, 1:im_h);
            M(:,x2-d) = 1;
        end

        % Find closest known point for LEFT boundary ------------------------------

        if settings.display == 1
            waitbar(0.675, wh, 'Decoding: Predicting borders');
        end

        ci = 1;
        Mtemp = M;
        Mtemp(:,ci+1:end) = 0;

        while sum(Mtemp(:)) < im_w/4
            ci = ci + 1;
            Mtemp = M;
            Mtemp(:,ci+1:end) = 0;
        end

        [known_y, known_x] = find(Mtemp ~= 0);
        x2 = 1;

        % imagesc(S.*Mtemp); hold on; xlim([1 ci]);
        for y2 = y_sam
            % find closest known point for bottom boundary
            best_distance = Inf;
            pixel_value = NaN;
            for i = 1:length(known_x)
                y1 = known_y(i);
                x1 = known_x(i);
                distance = sqrt((x1-x2)^2 + (y1-y2)^2);
                if distance < best_distance
                    best_distance = distance;
                    pixel_value = S(y1, x1);
        %             x1m = x1;
        %             y1m = y1;
                end
            end
            S(y2,x2) = pixel_value;
        end

        samples = S(y_sam, x2);
        for d = 0:border_dmax
            S(:,x2+d) = spline(y_sam, samples, 1:im_h);
            M(:,x2+d) = 1;
        end
        
    end

    if settings.display == 1
        waitbar(0.7, wh, 'Decoding: Pixel fusion completed');
    end

        
    % ----------------- Begin Iterative Filtering -----------------------------

    ker = @(a) [0.25-a/2, 0.25, a, 0.25, 0.25-a/2];
    krn = ker(0.4);
    Kernel = krn'*krn;

    Y = ones(size(S));

    L = settings.filter_levels;

    while L >= 1
        iters = 0;
        
        Y_old = zeros(size(S));
        
        fprintf('Level %d: ', L);
        if settings.display == 1
            waitbar(0.7 + (settings.filter_levels - L + 1)/settings.filter_levels*0.3, wh, 'Decoding: iterative filtering');
        end

        while mean2(abs(Y - Y_old)) > 1e-4
            
            
            if mod(iters, 20) == 0
                fprintf('(%.2f)', log10(mean2(abs(Y - Y_old))));
            end
            
            fprintf('.');

            Y_old = Y;
            Y = Y_old.*(1-M) + S.*M;

            for l = 1:L
                Y = conv2(krn, krn', Y, 'same');
                Y = Y(1:2:end, 1:2:end);        
            end

            for l = 1:L
                Yt = zeros(2*size(Y));
                Yt(1:2:end, 1:2:end) = Y;
                Y = 4*conv2(krn, krn', Yt, 'same');
            end

            iters = iters + 1;

        end
        
        if mod(iters-1, 20) ~= 0
            fprintf('(%.2f)', log10(mean2(abs(Y - Y_old))));
        end    
        fprintf(' - %d iters\n', iters);
        L = L - 1;
    end

    if settings.display == 1
        waitbar(1.0, wh, 'Decoding: decoding completed');
        close(wh);
    end


    Y = im2uint8(Y);
    error_map = reshape(error_map, size(Y)/32);

end
