function [outStruct] = setParamsWienerTest(structName,varargin)
%this function creates the input structures needed for performing
% the Wiener filter computation

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

switch structName
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

