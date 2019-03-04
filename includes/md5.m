function h = md5(inp)

  inp = inp(:);

  if ischar(inp) || islogical(inp)
      inp = uint8(inp);
  else
      inp = typecast(inp,'uint8');
  end

  % Get Java Message Digest
  md = java.security.MessageDigest.getInstance('MD5');
  md.update(inp);
  h = typecast(md.digest,'uint8');
  h = dec2hex(h)';
  h = lower(h(:)');
  clear x;
end