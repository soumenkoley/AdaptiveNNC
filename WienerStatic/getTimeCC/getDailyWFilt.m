% this script was written to observe all the cross-correlations and
% average them
% Next the matrices are used to compute the daily Wiener filter

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

clear; close all;

fPath = '/users/koley/CEBNNAnalysis/WienAdaptTilt/WienerStatic/getTimeCC/Output/';

gpsStart = 1374796818; % July 31, 2023
winLen = 1000;
nWin = 86;

RxxAvg = 0;
PxyAvg = 0;

for i = 1:1:nWin
    
    tNow = gpsStart + (i-1)*winLen;
    fName = [num2str(tNow),'.mat'];
    load([fPath,fName]);
    
    RxxAvg = RxxAvg+Rxx;
    PxyAvg = PxyAvg+Pxy;
    
%     figure(1);
%     hold on;
%     plot(Rxx(:,140));
%     hold off;
    
%     figure(2);
%     hold on;
%     plot(Pxy(:,13));
%     hold off;
end

RxxAvg = RxxAvg/nWin;
PxyAvg = PxyAvg/nWin;

figure(1);
hold on;
plot(RxxAvg(:,140),'k','LineWidth',2);
hold off;

figure(2);
hold on;
plot(PxyAvg(:,13),'k','LineWidth',2);
hold off;

% create the Toeplitz matrix
nRefCha = 24;
firOrder = 100;
R = zeros(nRefCha*(firOrder+1), nRefCha*(firOrder+1));

ccCount = 1;
for s1 = 1:1:nRefCha
    disp(['Estimating correlations for reference sensor ',num2str(s1)]);
    for s2 = s1:nRefCha
        rx = RxxAvg(:,ccCount);
        Rx = toeplitz(flipud(rx(1:(firOrder+1),1)), rx((firOrder+1):(2*firOrder+1),1));
        
        top = (s1-1)*(firOrder+1)+1;
        bottom = s1* (firOrder+1);
        left = (s2-1)*(firOrder+1)+1;
        right = s2*(firOrder + 1);
        
        R(top:bottom, left:right) = Rx;
        ccCount = ccCount+1;
        if(s1~=s2)
            R(left:right, top:bottom) = Rx'; % transpose
        end
    end
end

P = zeros(nRefCha*(firOrder+1),1);

ccCount = 1;
for s = 1:nRefCha
    px = PxyAvg(:,ccCount);
    
    top = (s-1)*(firOrder + 1) + 1;
    bottom = s*(firOrder+1);
    
    P(top:bottom,1) = px((firOrder+1):(2*firOrder+1),1);
    
    ccCount = ccCount+1;
end
rankOut = rank(R');
h = lsqminnorm(R',P);
%[Q,MM] = qr(R',0); QR decomposition will also give the same result
%gg = Q*(MM'\P);
hh = reshape(h,firOrder+1,nRefCha);
disp('done');