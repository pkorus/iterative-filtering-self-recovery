function crc = crc9(msg)
    % crc = MD5(uint8(msg), uint8(1));
    crc = md5(uint8(msg));
    crc = dec2binvec(mod(sum(uint8(crc)), 255), 8);
    
end
