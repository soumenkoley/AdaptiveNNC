function [chaName,nChannels] = readChannelNames(fName)
% read the station names to be read

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

fId = fopen(fName);

% this is the first line, no line skip needed
tLine = fgetl(fId);

tabInd = strfind(tLine, sprintf('\t'));

chaNo = 1;
chaName{chaNo} = tLine(1:(tabInd-1));

while(~feof(fId))
    chaNo = chaNo+1;
    tLine = fgetl(fId);
    tabInd = strfind(tLine, sprintf('\t'));
    chaName{chaNo} = tLine(1:(tabInd-1));
end

fclose(fId);

nChannels = length(chaName);
end