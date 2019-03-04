function [Xq, ci] = quantize_scalar(X, codebook)

      ci = uint16(1); si = uint16(1); fi = uint16(length(codebook));
      while fi - si > 1
        ci = bitshift(fi+si,-1);
        if X - codebook(ci) < 0
            fi = ci;
        else
            si = ci;
        end
      end
      if abs(X-codebook(si)) > abs(X-codebook(fi))
          ci = fi;
      else
          ci = si;
      end
      Xq = codebook(ci);

end