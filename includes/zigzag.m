function indices = zigzag(block)

% indices = zigzag(block)
%
% Returns indices of AC coefficients for watermark storage in a zig-zag
% manner. In case of block operation mode, the indices are selected in a
% per one per block manner. Eg. for block = 8
%
%   1 2 4 7 . . . . 
%   3 5 8 . . . . .
%   6 9 . . . . . .
%   . . . . . . . . 
%   . . . . . . . .
%   . . . . . . . .
%   . . . . . . . .
%   . . . . . . . .

indices = zeros(1,block^2);

i = 1;
for sc = 1:2*block
    cc = sc;
    cr = 1;
    while cc > 0
        if cc <= block && cr <= block
            indices(i) = (cc-1)*block+cr;
            i = i+1;
        end
        cc = cc - 1;
        cr = cr + 1;
    end
end

end