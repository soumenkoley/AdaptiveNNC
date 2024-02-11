function [yDec] = firDecimate(x,decFact,b,padYes)
% this function first applies a low pass fir filter to the data vector 'x'
% by convolving 'x' with 'b' (numerator coefficients of fir filter)
% padLen is the number of samples to be padded before
% and after the input data vector 'x'
% padYes = 0 or 1, 1 implies pad data otherwise 0
% at the last stage the decimation is performed using the interger
% factor 'decFact'

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

xLen = length(x);
newLen = floor(xLen/decFact);
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

yDec = zeros(newLen,1);

for m = 0:(newLen-1)
    yDec(m+1,1) = ySync(m*decFact+1,1);
end

end