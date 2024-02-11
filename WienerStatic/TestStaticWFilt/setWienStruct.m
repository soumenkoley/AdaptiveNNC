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
wStruct.compFlag = 1;
end

