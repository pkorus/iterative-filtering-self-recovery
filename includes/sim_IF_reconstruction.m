function [Y, psnr, psnr_m, ssim, S, line_distances] = sim_IF_reconstruction(I, line_bits, lines_count, prng_seed, settings)
% 
% [Y, psnr, psnr_m, ssim] = sim_IF_reconstruction(I, line_bits, lines_count, prng_seed, settings)
%
% Simulates iterative-filter-based reconstruction with a given number of
% lines of given precision.
%
% Params:
%
%   I                  - input image, double prec. numbers in range (0,1)
%
%   line_bits          - line representation precision (8 bits will be
%                        subtracted for descriptor hash)
%
%   lines_count        - number of lines to be sampled from the image
%
%   prng_seed          - seed for the pseudo random number generator
%
%   settings           - struct with advanced settings:
%     bit_alloc_metric - 1 (L1 distance), 2 (L2 distance)
%     line_selection   - <to be defined>
%     fusion_type      - 0 (progressive), 1 (best), 2 (mean), 3 (weighted)
%     disable_borders  - set 1 to disable border predictions
%     filter_levels    - levels of decomposition for iterative filtering
%
%
    
    if nargin < 5
        settings = default_settings();
    end

    % Generate random samples

    addpath('includes');
    if ~exist('codebook', 'var')
        load codebooks.mat
    end
    
    if settings.display == 1
        wh = waitbar(0.0, 'Simulation: Initialization');
    end

    [h, w] = size(I);

    ker = @(a) [0.25-a/2, 0.25, a, 0.25, 0.25-a/2];
    krn = ker(0.4);
    Kernel = krn'*krn;

    rand('twister', prng_seed);

    % Define representation fidelity for successive coefficients
    % q_map = [4 3*ones(1,10) 2*ones(1,6) 1*ones(1,6)];
    
    HASH_LEN = 8;
    if settings.bit_alloc_metric == 2
        q_map = allocate_line_descriptor_bits_L2(line_bits - HASH_LEN);
    else
        q_map = allocate_line_descriptor_bits(line_bits - HASH_LEN);
    end

    S = zeros(size(I));
    M = zeros(size(I));
    
    if settings.fusion_type > 1
        pixel_candidates = cell(size(I));
        pixel_weights = cell(size(I));
%     pixel_sources = cell(size(I));
    end

    line_distances = zeros(1,numel(lines_count));
    for k = 1:lines_count

        best_sum = Inf;
    %     best_sum = 0;
        coords_best = zeros(4,1);
        
        if settings.display == 1
            waitbar(k/lines_count*0.75, wh, 'Simulation: Processing successive descriptors');
        end

        if settings.profile_mode == 0
            % Plain, random selection of (x,y) -- (x,y)
            coords = rand(4,1);
            x1 = max(1, round(coords(1)*w));
            x2 = max(1, round(coords(2)*w));
            y1 = max(1, round(coords(3)*h));
            y2 = max(1, round(coords(4)*h));
            
            while (y2 == y1 && x2 == x1)
                coords = rand(4,1);
                x1 = max(1, round(coords(1)*w));
                x2 = max(1, round(coords(2)*w));
                y1 = max(1, round(coords(3)*h));
                y2 = max(1, round(coords(4)*h));
            end
                        
        elseif settings.profile_mode < 0
            % Select (x,y) then degree, then distance (according to
            % distrib)
            coords = rand(4,1);
            x1 = max(1, round(coords(1)*(w-10)));
            y1 = max(1, round(coords(2)*(h-10)));
            
            deg = linspace(0,2*pi,60);
            deg = deg(1:end-1);

            distances = zeros(1,length(deg));
            for d = deg

                [x3, y3] = find_coords(x1, y1, d, size(S,1), size(S,2));
            %     plot(x2, y2, 'go', 'MarkerSize', 10);
            %     title(sprintf('%d, %d', int32(x2), int32(y2)));
                distances(deg == d) = sqrt( (x1 - x3)^2 + (y1 - y3)^2 );
            end

            % plot(deg, distances);

            min_d = min(distances);
            max_d = max(distances);

            % Adjust this line to control the distribution of lengths
            new_d = min_d + (0.05 + 0.95*coords(3)^5)*(max_d - min_d);
            valid_degs = find(distances >= new_d);
            final_deg = deg(valid_degs(1 + round(coords(4)*(numel(valid_degs)-1))));

            [x3, y3] = find_coords(x1, y1, final_deg, size(S,1), size(S,2));
            dx = x3 - x1;
            dy = y3 - y1;
            x2 = round(x1 + dx*new_d/sqrt( (x1 - x3)^2 + (y1 - y3)^2 ));
            y2 = round(y1 + dy*new_d/sqrt( (x1 - x3)^2 + (y1 - y3)^2 ));

        elseif settings.profile_mode >= 1
            % Best of N lines
            for l = 1:settings.profile_mode
                coords = rand(4,1);
                x1 = max(1, round(coords(1)*w));
                x2 = max(1, round(coords(2)*w));
                y1 = max(1, round(coords(3)*h));
                y2 = max(1, round(coords(4)*h));

                while (y2 == y1 && x2 == x1)
                    coords = rand(4,1);
                    x1 = max(1, round(coords(1)*w));
                    x2 = max(1, round(coords(2)*w));
                    y1 = max(1, round(coords(3)*h));
                    y2 = max(1, round(coords(4)*h));
                end

                new_sum = sqrt((x1-x2)^2 + (y1-y2)^2);
                if new_sum < best_sum
                    coords_best = coords;
                    best_sum = new_sum;
                end
            end
            
            x1 = max(1, round(coords_best(1)*w));
            x2 = max(1, round(coords_best(2)*w));
            y1 = max(1, round(coords_best(3)*h));
            y2 = max(1, round(coords_best(4)*h));
            
        end

        line_distances(k) = sqrt((x1-x2)^2 + (y1-y2)^2);

        x = linspace(x1,x2,1000);
        y = linspace(y1,y2,1000);
        index = sub2ind(size(I),round(y),round(x));    
        scan_profile = I(index);
        scan_spectrum = dct(scan_profile);


        % quantize coefficients
        % 1. dc 7 bits uniform
        % 2. ac 6/5/4 bits Lloyd
        scan_spectrumq = zeros(size(scan_spectrum));
        for i = 1:4
            [scan_spectrumq(q_map == i), ~] = quantize(scan_spectrum(q_map == i), codebook{i});
        end
        
        scan_profile_compq = idct(scan_spectrumq);
        
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
        
    if settings.fusion_type > 1
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

    if settings.disable_borders ~= 1
    
        % Try to predict borders
        if settings.display == 1
            waitbar(0.75, wh, 'Simulation: Predicting borders');
        end

        border_dmax = 3;

        [im_h, im_w] = size(S);

        % Sample each border using 32 points - put value from closes known point
        x_sam = round(linspace(1, im_w, 64));
        y_sam = round(linspace(1, im_h, 64));

        % Find closest known point for UPPER boundary

        ci = 1;
        Mtemp = M;
        Mtemp(ci+1:end,:) = 0;
        known_indices = find(Mtemp ~= 0);
        while length(known_indices) < im_h
            ci = ci + 1;
            Mtemp = M;
            Mtemp(ci+1:end,:) = 0;
            known_indices = find(Mtemp ~= 0);
        end

        known_indices = known_indices';

        y2 = 1;
        for x2 = x_sam
            best_distance = Inf;
            pixel_value = NaN;
            for i = known_indices
                [y1, x1] = ind2subq([im_h im_w], i);
                distance = sqrt((x1-x2)^2 + (y1-y2)^2);
                if distance < best_distance
                    best_distance = distance;
                    pixel_value = S(i);
                end
        %         delete(h);
            end
            S(y2,x2) = pixel_value;
        end

        samples = S(y2, x_sam);
        for d = 0:border_dmax
            S(y2+d,:) = spline(x_sam, samples, 1:im_w);
            M(y2+d,:) = 1;
        end

        % Find closest known point for BOTTOM boundary
        if settings.display == 1
            waitbar(0.7625, wh, 'Simulation: Predicting borders');
        end

        ci = 1;
        Mtemp = M;
        Mtemp(1:end-ci,:) = 0;
        known_indices = find(Mtemp ~= 0);
        while length(known_indices) < im_h
            ci = ci + 1;
            Mtemp = M;
            Mtemp(1:end-ci,:) = 0;
            known_indices = find(Mtemp ~= 0);
        end

        known_indices = known_indices';

        y2 = im_h;

        for x2 = x_sam
            % find closest known point for bottom boundary
            best_distance = Inf;
            pixel_value = NaN;
            for i = known_indices
                [y1, x1] = ind2subq([im_h im_w], i);
                distance = sqrt((x1-x2)^2 + (y1-y2)^2);
                if distance < best_distance
                    best_distance = distance;
                    pixel_value = S(i);
                end
            end
            S(y2,x2) = pixel_value;
        end

        samples = S(y2, x_sam);
        for d = 0:border_dmax
            S(y2-d,:) = spline(x_sam, samples, 1:im_w);
            M(y2-d,:) = 1;
        end

        % Find closest known point for RIGHT boundary
        if settings.display == 1
            waitbar(0.775, wh, 'Simulation: Predicting borders');
        end

        ci = 1;
        Mtemp = M;
        Mtemp(:,1:end-ci) = 0;
        known_indices = find(Mtemp ~= 0);
        while length(known_indices) < im_w
            ci = ci + 1;
            Mtemp = M;
            Mtemp(:,1:end-ci) = 0;
            known_indices = find(Mtemp ~= 0);
        end

        known_indices = known_indices';

        x2 = im_w;

        for y2 = y_sam
            % find closest known point for bottom boundary
            best_distance = Inf;
            pixel_value = NaN;
            for i = known_indices
                [y1, x1] = ind2subq([im_h im_w], i);
                distance = sqrt((x1-x2)^2 + (y1-y2)^2);
                if distance < best_distance
                    best_distance = distance;
                    pixel_value = S(i);
                end
            end
            S(y2,x2) = pixel_value;
        end

        samples = S(y_sam, x2);
        for d = 0:border_dmax
            S(:,x2-d) = spline(y_sam, samples, 1:im_h);
            M(:,x2-d) = 1;
        end

        % Find closest known point for LEFT boundary
        if settings.display == 1
            waitbar(0.7825, wh, 'Simulation: Predicting borders');
        end
        
        ci = 1;
        Mtemp = M;
        Mtemp(:,ci+1:end) = 0;
        known_indices = find(Mtemp ~= 0);
        while length(known_indices) < im_w
            ci = ci + 1;
            Mtemp = M;
            Mtemp(:,ci+1:end) = 0;
            known_indices = find(Mtemp ~= 0);
        end

        known_indices = known_indices';

        x2 = 1;

        for y2 = y_sam
            % find closest known point for bottom boundary
            best_distance = Inf;
            pixel_value = NaN;
            for i = known_indices
                [y1, x1] = ind2subq([im_h im_w], i);
                distance = sqrt((x1-x2)^2 + (y1-y2)^2);
                if distance < best_distance
                    best_distance = distance;
                    pixel_value = S(i);
                    x1m = x1;
                    y1m = y1;
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
    
    % Perform reconstruction
    Y = ones(size(S));

    L = settings.filter_levels;
    
    if settings.display == 1
            waitbar(0.80, wh, 'Simulation: Iterative filtering');
    end

    while L >= 1
        iters = 0;
        Y_old = zeros(size(S));
        
        if settings.display == 1
            waitbar(0.80 + 0.2*(settings.filter_levels-L)/settings.filter_levels, wh, 'Simulation: Iterative filtering');
        end

%         fprintf('Level %d: ', L);

        while mean2(abs(Y - Y_old)) > 1e-4

            Y_old = Y;
            Y = Y_old.*(1-M) + S.*M;

            for l = 1:L
                Y = conv2(Y,Kernel, 'same');
                Y = Y(1:2:end, 1:2:end);        
            end

            for l = 1:L
                Yt = zeros(2*size(Y));
                Yt(1:2:end, 1:2:end) = Y;
                Y = 4*conv2(Yt, Kernel, 'same');
            end

            iters = iters + 1;

        end

        L = L - 1;
    end
    
    if settings.display == 1
        waitbar(1.0, wh, 'Simulation: Completed');
        close(wh);
    end

    mask = ones(size(Y));
    mask(1:4,:) = 0;
    mask(end-4:end,:) = 0;
    mask(:,1:4) = 0;
    mask(:,end-4:end) = 0;

    psnr_m = calc_psnr_map(im2uint8(I), im2uint8(Y), mask);
    psnr = calc_psnr(im2uint8(I), im2uint8(Y));
    ssim = ssim_index(im2uint8(I), im2uint8(Y));
%     fprintf('PSNR: %.2f\n', psnr);

    imwrite(M, 'last_mask.png');
    imwrite(S, 'last_input.png');
    imwrite(Y, 'last_result.png');
end
