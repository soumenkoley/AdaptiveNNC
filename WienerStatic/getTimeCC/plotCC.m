% this script was written to observe all the cross-correlations and
% how they change in magnitude per 1000 s and after avergaing over a day

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
    
    figure(1);
    hold on;
    plot(Rxx(:,140));
    hold off;
    
    figure(2);
    hold on;
    plot(Pxy(:,13));
    hold off;
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