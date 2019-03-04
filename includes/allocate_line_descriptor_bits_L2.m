function [q_map, dist_L1, dist_L2] = allocate_line_descriptor_bits_L2(bits)

%     if nargin < 2
%         mode = 1;
%     end
%     
%     if mode == 1
%         load dists.mat
%     elseif mode == 2
%         load dists_L1_new.mat
%     elseif mode == 3
%         load dists_L2_new.mat
%     end
    
    load dists.mat
    
    q_map = zeros(1,1000);
    q_map(1) = 4;

    allocated = 6*sum(q_map==4) + 6*sum(q_map==3) + 5*sum(q_map==2) + 4*sum(q_map==1);

    while allocated < bits
        % Find location, where allocation upgrade will benefit the most
        base_d_map = distances(sub2ind(size(distances), q_map(2:end)+1, 1:999));
        base_dist = sqrt(sum(base_d_map.^2));
        update_map = zeros(1,1000);
        bit_change = ones(1,1000);
        for i = 2:1000
            q_mapu = q_map;
            q_mapu(i) = q_mapu(i) + 1;
            
            if (q_mapu(i) == 1 && bits - allocated < 4) || q_mapu(i) == 4
                continue
            end
            
            allocated = 6*sum(q_mapu==4) + 6*sum(q_mapu==3) + 5*sum(q_mapu==2) + 4*sum(q_mapu==1);
            
            if (bits - allocated < 4 && bits - allocated > 2*sum(q_mapu==1) + sum(q_mapu==2))
                continue
            end
            
%             if (q_mapu(i) == 1 && bits - allocated < 4) || q_mapu(i) == 4
%                 continue
%             end

%             new_d_map = base_d_map;
%             new_d_map(i) = distances( % distances(sub2ind(size(distances), q_mapu(2:end)+1, 1:999));
            
%             update_map(i) = abs(base_d_map(i-1) - distances(q_mapu(i)+1, i-1));
            
%             update_map(i) = sum(abs(base_d_map - distances(sub2ind(size(distances), q_mapu(2:end)+1, 1:999))));
%             update_map(i) = sqrt(sum(abs(base_d_map - distances(sub2ind(size(distances), q_mapu(2:end)+1, 1:999))).^2));

            update_map(i) = base_dist - sqrt(sum(abs(distances(sub2ind(size(distances), q_mapu(2:end)+1, 1:999))).^2));
            
%             update_map(i) = sum(base_d_map) - sum(abs(distances(sub2ind(size(distances), q_mapu(2:end)+1, 1:999))));
            
            if q_mapu(i) == 1
                bit_change(i) = 4;
            end
        end
        % Find best change place
        change_map = update_map./bit_change;
        [~, indices] = sort(change_map);
        % Make the change
        q_map(indices(end)) = q_map(indices(end)) + 1;
        allocated = 6*sum(q_map==4) + 6*sum(q_map==3) + 5*sum(q_map==2) + 4*sum(q_map==1);
    end

    if allocated ~= bits
        fprintf('ERROR Allocated %d bits out of %d \n', allocated, bits);
    end

    dist_L2 = sqrt(sum(distances(sub2ind(size(distances), q_map(2:end)+1, 1:999)).^2));
    dist_L1 = sum(distances(sub2ind(size(distances), q_map(2:end)+1, 1:999)));
    
    maxi = find(q_map, 1, 'last');
    q_map = q_map(1:maxi);

end