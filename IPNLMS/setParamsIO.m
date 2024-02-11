function [outStruct] = setParamsIO(structName,varargin)
%this function creates the input structures needed for running
% several of the main codes

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

switch structName
    case 'readStruct'
        % the inputs for reading the data
        outStruct.gpsStart = 1374883218; % Aug 01, 2023
        % make sure the gps time is a multiple of 100
        % gwf files are stored every 100 s
        outStruct.sigLen = 1000; % read sigLen in one go, unit in secs
        outStruct.nWin = 86*7; % number of 1000 s windows
        % buffer is necessary because of the application of the FIR
        % filters that introduces a delay
        outStruct.buffLen = 5; % unit is seconds
        outStruct.removeLen = 4;
        % bufflen is added to the start and end of signal
        % to reduce the effects of filtering on the data
        outStruct.saveData = 1; % if 1, then save the data, else 0
        outStruct.savePath = '0'; % some dummy string
        outStruct.saveName = [num2str(outStruct.gpsStart)];
        outStruct.nFiles = floor(outStruct.sigLen/100)-1;
        outStruct.gwfFileName = 'gwfFName.txt'; % this is the file
        % where all the gwf file paths are stored
    case 'procStruct'
        % reads the inputs used for pre-processing the data
        outStruct.fSampRef = 500; % sampling frequency for reference channels
        outStruct.fSampTar = 10000; % sampling frequecy for target channel
        
        % low pass filter before decimation
        outStruct.lpOrdRef = 51; % order of the lowpass filter for the geophones
        outStruct.lpOrdTar = 5001; % order of the lowpass filter for the tilt-signal
        
        outStruct.fLowCut = 30; % units in Hz of low-pass frequency
        
        % high-pass filter after decimation
        outStruct.hpOrd = 51; % order for the high pass filter
        outStruct.fHighCut = 10; % units in Hz of the high-pass frequency

        % downsampling frequency
        outStruct.fDownSamp = 100; % units in Hz

        % decimation factor
        outStruct.dsFactRef = floor(outStruct.fSampRef/outStruct.fDownSamp);
        outStruct.dsFactTar = floor(outStruct.fSampTar/outStruct.fDownSamp);
        % flag for padding to be done before filtering
        outStruct.padYes = 1;

        % bias of the reference and target channels
        outStruct.tarBias = 6.107440109204276e-10;
        outStruct.refBias = 4.869999958856397e-10;

        % flag to turn on scaling of data
        outStruct.scaleInp = 0; % default to 1
    case 'outDataStruct'
        nRefCha = varargin{1};
        readStruct = varargin{2};
        procStruct = varargin{3};
        nSamp = (readStruct.sigLen + 2*readStruct.buffLen)*procStruct.fDownSamp;
        outStruct.refData  = zeros(nSamp,nRefCha);
        outStruct.tarData = zeros(nSamp,1);
        outStruct.fSamp = procStruct.fDownSamp;
        outStruct.saveName = '0';
        outStruct.ifPlot = 0; % set to 1, if you want to see the 
        % wiener output and the error
        outStruct.scaleInp = procStruct.scaleInp; % default to 1
     
    otherwise
        error('Incorrect input structure name');
end

end

