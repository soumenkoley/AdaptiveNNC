function [frls] = setFRLSStruct(firOrder,nChannels)
% This script sets up the data structure for the FRLS algorithm

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

P = nChannels;
L = firOrder;
PL = P*L;

frls.start = 1;
frls.firstEntry = 1;
frls.lambda = 1 - 1/(3*PL);
frls.lambdaN = frls.lambda^(PL);
frls.lambda1 = 1/frls.lambda;
frls.mu = 5*10^10;
frls.gamma1 = 1.0;
frls.gamma = 1.0;
% feedback coefficients for numerical stabilization
frls.K = [1.5,2.5,1,0,1,1];
frls.fSamp = 100; % data was downsampled to 100 Hz
frls.sigLen = 1000; % units in seconds
frls.nSamp = frls.fSamp*frls.sigLen;

% forward predictor
frls.A = [eye(P),zeros(P,PL)];
% backward predictor
frls.B = [zeros(P,PL),eye(P)];

% kalman gain
frls.C = zeros(PL,1);
frls.CN1 = zeros((PL+P),1);

% filter coefficients
frls.W = zeros(PL,1);

% reference channel data
frls.X = zeros((PL+P),1);
frls.XX = zeros((PL+P),1);

% other matrix partitions 
frls.alpha1 = (1/(frls.lambdaN*frls.mu))*eye(P);
frls.beta = frls.mu*eye(P);
frls.betaInv = 1/frls.mu*eye(P);

% rescue must always be turned on at start
frls.doRescue = 1;
end
