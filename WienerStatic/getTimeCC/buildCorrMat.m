% this script was written to build the daily correlation matrices
% that will be used to compute the correlation matrix
% updated Aug 24, 2023

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

clear; %close all;

% addpath to several subroutines
run  /users/swinkels/deploy/MatlabVirgoTools/trunk/startup.m

currDir = pwd;

firOrder = 101;

% initial the read Data structure
[outReadStruct] = setParamsIO('readStruct');

outReadStruct.savePathWiener = [currDir,'/Output/'];

% read the reference and target channel names
[chaName,nChannels] = readChannelNames('RefChaNames.txt');

% create the lowpass and the highpass filters for preprocessing
% initialize the pre-processing data structure
procStruct = setParamsIO('procStruct');
procStruct.scaleInp = 0; % turn this on while scaling

% prepare the Wiener input structure
wienStruct = setWienStruct(firOrder);

% get the numerator coefficients for the lowpass filters
% they are created separately for target and reference since
% they are sampled at different frequencies
[hLPTar] = getHammLPFir(procStruct.lpOrdTar,procStruct.fLowCut,...
     procStruct.fSampTar);
% 
[hLPRef] = getHammLPFir(procStruct.lpOrdRef,procStruct.fLowCut,...
     procStruct.fSampRef);
% 
% % the highpass filters are same for both target and reference
[hHP] = getHammHPFir(procStruct.hpOrd,procStruct.fHighCut,...
     procStruct.fDownSamp);
% 
% % lowpass again with reduced sampling rate
[hLPNew] = getHammLPFir(procStruct.lpOrdRef,procStruct.fLowCut,...
     procStruct.fDownSamp);

% loop through all channel names, pre-process data and store
% for later use

% check if the correlation matrices are full rank
rankMat = zeros(outReadStruct.nWin,1);

for winNo = 1:1:outReadStruct.nWin
    % initialize the output data structure
    outDataStruct = setParamsWienerTest('outDataStruct',nChannels,outReadStruct,...
                                   procStruct);
    startTimeStr = num2str(outReadStruct.gpsStart + ...
        outReadStruct.sigLen*(winNo-1));
    
    outDataStruct.saveNameWiener = [outReadStruct.savePathWiener,startTimeStr,'.mat'];
    
    startTime = outReadStruct.gpsStart + ...
        outReadStruct.sigLen*(winNo-1)-outReadStruct.buffLen; % 2secs buffer, get rid
    % of filtering effects at the edges
    sigLen = outReadStruct.sigLen + 2*outReadStruct.buffLen; % buffer at the end too
    
    % read the target channel first
    yTemp = getTiltSig(startTime,sigLen,0,0);
    % scale by precomputed scaking factor
    yTemp = yTemp/procStruct.tarBias;
    [yDec] = firDecimate(yTemp,procStruct.dsFactTar,hLPTar,...
                 procStruct.padYes);
    outDataStruct.tarData = yDec;
    yTemp =[]; yDec =[];
    disp('Target data read!');
    
    for chaNo = 1:1:nChannels
       
        dataNow = getChannel(chaName{chaNo},startTime,sigLen);
        dataNow = dataNow/procStruct.refBias;
        if(chaNo>=1)
            % these are the reference channels
            
            [yDec] = firDecimate(dataNow,procStruct.dsFactRef,hLPRef,...
            procStruct.padYes);
            
            if(procStruct.scaleInp==1)
                % remove the mean
                yDec = yDec-mean(yDec);
            end
            
            outDataStruct.refData(:,chaNo) = yDec;
            
        else
           disp('Error loading reference channels');
        end
    end

    % rescale the data by standard deviation of the target and
    % reference data
    if(winNo==1)
        stdRef = std(outDataStruct.refData(:));
        stdTar = std(outDataStruct.tarData);
    else
        if(wienStruct.compFlag==1)
            stdRef = std(outDataStruct.refData(:));
            stdTar = std(outDataStruct.tarData);
        end
    end
    
    if(procStruct.scaleInp==1)
        % rescale data
        outDataStruct.refData = outDataStruct.refData/stdRef;
        outDataStruct.tarData = outDataStruct.tarData/stdTar;
    end

    % apply the high pass and low pass again
    for i = 1:1:length(outDataStruct.refData(1,:))
        outDataStruct.refData(:,i) = firFilt(outDataStruct.refData(:,i),hHP,...
            procStruct.padYes);
        outDataStruct.refData(:,i) = firFilt(outDataStruct.refData(:,i),hLPNew,...
            procStruct.padYes);
    end

    outDataStruct.tarData = firFilt(outDataStruct.tarData,hHP,...
        procStruct.padYes);
    outDataStruct.tarData = firFilt(outDataStruct.tarData,hLPNew,...
        procStruct.padYes);

    % the buffer was 5 s on either end of the signal
    % we remove 4 s in start and 4 s at the end 
    % and compute the Wiener filter
    % initially the signal was from -5 to 1005,
    % now its from -1 to 1001
    % removing one second at the start and end is enough
    % to care of edge effects
    lenInit = length(outDataStruct.tarData);
    startBuff = outReadStruct.removeLen*procStruct.fDownSamp+1;
    endBuff = lenInit - outReadStruct.removeLen*procStruct.fDownSamp;
    outDataStruct.refData = outDataStruct.refData(startBuff:endBuff,:);
    outDataStruct.tarData = outDataStruct.tarData(startBuff:endBuff,1);
    
    % build the correlation matrix
    if(winNo==1)
        % first time entry, compute Wiener filter
        [Rxx,Pxy] = wienerMISO(outDataStruct,wienStruct.firOrder-1);
    else
        if(wienStruct.compFlag==1)
            [Rxx,Pxy] = wienerMISO(outDataStruct,wienStruct.firOrder-1);
        end
    end

    save(outDataStruct.saveNameWiener ,'Rxx','Pxy');

    disp('1000 s done!');

end

% end of code
