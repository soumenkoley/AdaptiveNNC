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
% this must be set to 1.0 for recomputing the Wiener filter
% every 1000 s
wStruct.compFlag = 1;
end

