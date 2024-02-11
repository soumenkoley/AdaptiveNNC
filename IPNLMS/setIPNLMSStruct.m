function [ipnlms] = setIPNLMSStruct(P,L)
%this functions sets up the data structure for the
% IPNLMS algorithm
% update on Aug 22, 2023
%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

ipnlms.alpha = 1;
ipnlms.cst = 0.01;
ipnlms.ep = 10^-10;
ipnlms.deltaNLMS = 10^-5;
ipnlms.K = 0;
ipnlms.X = zeros(L,P);
ipnlms.H = zeros(L,P);
ipnlms.sigmaX = 0;
ipnlms.sigLen = 1000;
ipnlms.fSamp = 100;
ipnlms.nSamp = ipnlms.fSamp*ipnlms.sigLen;
ipnlms.sigmaBuff = 3; % units in seconds
ipnlms.startBuff = 1;
ipnlms.start = 1;

end