function [wStruct] = setWienStruct(firOrder)
%this function creates the data structure for the Wiener
% filter

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

wStruct.firOrder = firOrder;
% compFlag must be turned to 1, since we compare the dynamic Wiener
% performance with the FRLS algorithm
% essentially filter coefficients are recomputed every 1000 s
wStruct.compFlag = 1;
end

