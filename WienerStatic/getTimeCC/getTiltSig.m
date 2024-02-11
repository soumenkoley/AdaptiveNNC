function [F] = getTiltSig(gpsStart,sigLen,ifPlot,ifPlotTs)
% this function was written to get the calibrated output signal from
% the tiltmeter at NEB

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

% load the NNC_NEB_TILT_ITF
% A = getChannel('V1:NNC_NEB_TILT_ITF',gpsStart,sigLen);
B = getChannel('V1:NNC_NEB_TILT_PickOff',gpsStart,sigLen);

F = getChannel('V1:NNC_NEB_TILT_ERROR_PRE',gpsStart,sigLen);
F = F*3.8*10^-7; % do calibration

% multiply with pickoff
F = F.*B;
% now get the pwelch

if(ifPlot)
    [pxxF,fVec] = pwelch(F,200000,50000,200000,10000);
    %[pxxFDec,fVecDec] = pwelch(FDec,20*100,5*100,20*100,100);
    figure(1)
    hold on;
    plot(fVec,sqrt(pxxF),'b');
    %plot(fVecDec,sqrt(pxxFDec),'r')
    hold off;
end

if(ifPlotTs)
    figure(2);
    hold on;
    plot((1:1:length(F))/10000,F,'b');
    %plot((1:1:length(FDec))/100,FDec,'r')
    hold off;
end

%disp('Testing!');
end