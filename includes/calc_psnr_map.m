function psnr = calc_psnr_map(imga, imgb, map)
    imga = double(im2uint8(imga));
    imgb = double(im2uint8(imgb));
    psnr = 10*log10(65025/ (sum(sum(map.*((imga-imgb).^2)))/sum(map(:))));
end