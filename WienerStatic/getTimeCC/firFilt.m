function [ySync] = firFilt(x, b,padYes)
% this function applies a fir filter to the data vector 'x'
% but convolving 'x' with 'b' (numerator coefficients of fir filter)
% padLen is the number of samples to be padded before
% and after the input data vector 'x'
% padYes = 0 or 1, 1 implies pad data otherwise 0

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

xLen = length(x);
M = length(b);

if(padYes)
    padLen = 2*(M-1);
    y = padData(x,padLen);
else
    y = x;
end

y = conv(b,y);

if(padYes)
    stInd = padLen + floor((M-1)/2) + 1;
    ySync = y(stInd:(stInd+xLen-1),1);
else
    ySync = y;
end

end