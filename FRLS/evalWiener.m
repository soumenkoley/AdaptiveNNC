function [wienerOut,tarData] = evalWiener(inpStruct,wienerFilt)

% this function applies the Wiener filter to the data

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

[nSampFilt,nSens] = size(wienerFilt);
[nSampData,~] = size(inpStruct.refData);

wienerOutConv = zeros(nSampData+nSampFilt-1,1);

for i = 1:1:nSens
    wienerOutConv = wienerOutConv+conv(inpStruct.refData(:,i),wienerFilt(:,i));
end

wienerOut = wienerOutConv(nSampFilt:(length(wienerOutConv)-(nSampFilt-1)),1);
lWienerOut = length(wienerOut);
% now removing the last 100 samples
% remember, this output is for 1-1001 s
% so removing the last 1 second before comparing
wienerOut = wienerOut(1:(lWienerOut-nSampFilt+1),1);
% the target data as before was from -1 to 1001 s
% now remonving the first and the last 1 second
tarData = inpStruct.tarData(nSampFilt:(nSampData-nSampFilt+1),1);

wienErr = tarData-wienerOut;

if(inpStruct.ifPlot==1)
    figure(1);
    subplot(2,1,1);
    hold on;
    plot((1:1:length(tarData))/inpStruct.fSamp,tarData,'b');
    plot((1:1:length(tarData))/inpStruct.fSamp,wienerOut,'r');
    hold off;

    subplot(2,1,2);
    hold on;
    plot((1:1:length(tarData))/inpStruct.fSamp,wienErr,'r');
    hold off;
end

%save(inpStruct.saveName,'tarData','wienerOut');
disp('Wiener Applied');
end