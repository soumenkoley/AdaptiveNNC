function [y] = padData(x,padLen)
% this function pads the data with padLen samples on either side

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

xLen = length(x);
if(padLen>xLen)
    error('pad length cannot be greater than the input signal length');
end

y = zeros(xLen + 2*padLen,1);

y((padLen+1):(padLen+xLen),1) = x;
y(1:padLen,1) = flipud(x(1:padLen,1));

y((padLen+xLen+1):end,1) = flipud(x((xLen-padLen+1):end,1));

end