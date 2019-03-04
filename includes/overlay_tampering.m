function [ tamp_map ] = overlay_tampering(img, error_map)
    err = imresize(error_map, size(img), 'nearest');
    err1 = zeros(size(err));
    err1(err == 1) = 1;
    err1 = edge(err1);
    err1 = imdilate(err1, strel('disk',1));
    tamp_map = imoverlay(img, err1, [1 0 0]);
    err2 = zeros(size(err));
    err2(err == 2) = 1;
    err2 = edge(err2);
    err2 = imdilate(err2, strel('disk',1));
    tamp_map = imoverlay(tamp_map, err2, [0 0 1]);
end

