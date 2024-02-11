% this script was written to test the performance of the
% static Wiener filter compared to the dynamic Wiener filter

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

% load the Wiener filter calculated using a day of cross-correlations
load('WFiltJuly31.mat');
wienerFilt = hh;
hh=[];

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
    % scale the target data by precomputed scaling factor
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
            %yDec = decimate(dataNow,5);
            %yDec = filtfilt(HdBp10_30Hz.Numerator,1,yDec);
            if(procStruct.scaleInp==1)
                % remove the mean
                yDec = yDec-mean(yDec);
            end
            
            outDataStruct.refData(:,chaNo) = yDec;
            
        else
            
            disp('Error loading reference channels!');
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
    
    
    % apply the static wiener filter on the input channels
    [wienerOut,~] = evalWiener(outDataStruct,wienerFilt);
    
    if(winNo==1)
        % first time entry, compute Wiener filter
        [DWFilt,~] = wienerMISO(outDataStruct,wienStruct.firOrder-1);
    else
        if(wienStruct.compFlag==1)
            [DWFilt,~] = wienerMISO(outDataStruct,wienStruct.firOrder-1);
        end
    end

    % apply the dynamic wiener filter on the input channels
    [DWOut,tarData] = evalWiener(outDataStruct,DWFilt);
    
    % save the outputs
    %save(outDataStruct.saveNameWiener ,'wienerOut','tarData','wienerFilt');
    
    %do the ffts every 10 s and check, hard coded
    for mm = 1:1:100
        sm1 = (mm-1)*1000+1;
        sm2 = sm1+1000-1;
        [tFFT,f1] = pwelch(tarData(sm1:sm2,1),100,50,100,100);
        [wFFTErr,~] = pwelch(wienerOut(sm1:sm2,1)-tarData(sm1:sm2,1),100,50,100,100);
        [DWFFTErr,~] = pwelch(DWOut(sm1:sm2,1)-tarData(sm1:sm2,1),100,50,100,100);
        f1Ind = find(f1>=10,1,'first');
        f2Ind = find(f1>=15,1,'first');
        f3Ind = find(f1>=20,1,'first');
        f4Ind = find(f1>=25,1,'first');
        
        errAll10_15(mm,1) = 10*log10(mean(wFFTErr(f1Ind:f2Ind,1))/mean(tFFT(f1Ind:f2Ind,1)));
        errAll10_15(mm,2) = 10*log10(mean(DWFFTErr(f1Ind:f2Ind,1))/mean(tFFT(f1Ind:f2Ind,1)));
        
        errAll15_20(mm,1) = 10*log10(mean(wFFTErr(f2Ind:f3Ind,1))/mean(tFFT(f2Ind:f3Ind,1)));
        errAll15_20(mm,2) = 10*log10(mean(DWFFTErr(f2Ind:f3Ind,1))/mean(tFFT(f2Ind:f3Ind,1)));
        
        errAll20_25(mm,1) = 10*log10(mean(wFFTErr(f3Ind:f4Ind,1))/mean(tFFT(f3Ind:f4Ind,1)));
        errAll20_25(mm,2) = 10*log10(mean(DWFFTErr(f3Ind:f4Ind,1))/mean(tFFT(f3Ind:f4Ind,1)));
        
    end
    
    figure(1);
    subplot(1,3,1)
    plot((1:1:length(errAll10_15(:,1)))*10,errAll10_15(:,1),'b');
    hold on;
    plot((1:1:length(errAll10_15(:,1)))*10,errAll10_15(:,2),'r');
    ylim([-15,0]);
    ylabel('r_{10,15} (dB)');
    xlabel('Time (s)');
    title('Static vs DWF');
    legend({'Static Wiener','Dynamic Wiener'});
    hold off;
    
    subplot(1,3,2)
    plot((1:1:length(errAll15_20(:,1)))*10,errAll15_20(:,1),'b');
    hold on;
    plot((1:1:length(errAll15_20(:,1)))*10,errAll15_20(:,2),'r');
    ylim([-15,0]);
    ylabel('r_{15,20} (dB)');
    xlabel('Time (s)');
    title('Static vs DWF');
    legend({'Static Wiener','Dynamic Wiener'});
    hold off;
    
    subplot(1,3,3)
    plot((1:1:length(errAll20_25(:,1)))*10,errAll20_25(:,1),'b');
    hold on;
    plot((1:1:length(errAll20_25(:,1)))*10,errAll20_25(:,2),'r');
    ylim([-15,0]);
    ylabel('r_{20,25} (dB)');
    xlabel('Time (s)');
    title('Static vs DWF');
    legend({'Static Wiener','Dynamic'});
    hold off;
%     
    disp('1000 s done!');
    useInp = input('Do you wish to continue for the next 1000 s? Y or N \n','s');
    if(useInp=='N')
        break;
    end

end

% end of code