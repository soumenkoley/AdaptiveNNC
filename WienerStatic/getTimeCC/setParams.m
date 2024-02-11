function [outStruct] = setParams(structName,varargin)
%this function creates the input structures needed for running
% several of the main codes

switch structName
    case 'readStruct'
        % the inputs for reading the data
        outStruct.gpsStart = 1359417600; % gps time
        % make sure the gps time is a multiple of 100
        % gwf files are stored every 100 s
        outStruct.sigLen = 1000; % read sigLen in one go, unit in secs
        outStruct.saveData = 1; % if 1, then save the data, else 0
        outStruct.savePath = '0'; % some dummy string
        outStruct.saveName = [num2str(outStruct.gpsStart),'.mat'];
        outStruct.nFiles = floor(outStruct.sigLen/100)-1;
        outStruct.gwfFileName = 'gwfFName.txt'; % this is the file
        % where all the gwf file paths are stored
    case 'procStruct'
        % reads the inputs used for pre-processing the data
        outStruct.fSampRef = 500; % sampling frequency for reference channels
        outStruct.fSampTar = 1000; % sampling frequecy for target channel
        
        % low pass filter before decimation
        outStruct.lpOrd = 51; % order of the lowpass filter
        outStruct.fLowCut = 45; % units in Hz of low-pass frequency

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
    case 'outDataStruct'
        nRefCha = varargin{1};
        readStruct = varargin{2};
        procStruct = varargin{3};
        nSamp = readStruct.sigLen*procStruct.fDownSamp;
        outStruct.refData  = zeros(nSamp,nRefCha);
        outStruct.tarData = zeros(nSamp,1);

    otherwise
        error('Incorrect input structure name');
end

end

