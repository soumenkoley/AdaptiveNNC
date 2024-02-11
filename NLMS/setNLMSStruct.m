function [nlms] = setNLMSStruct(P,L)
%this functions sets up the data structure for the
% NLMS algorithm

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

nlms.firOrder = L;
nlms.alpha = 1;
nlms.delta = 10^-5;
nlms.sigLen = 1000;
nlms.filt = zeros(L,P);

end