function [nlmsStruct,xEst,err,filtOut] = doNLMS(inpStruct,nlmsStruct)
%this function iterates of samples of all the input channels
% and creates an adaptive filter following the steepest descent method
% nlmsOut is the final output
% nlmsError is the error between the measured and the estimated output

% inpStruct.tarData = inpStruct.tarData-mean(inpStruct.tarData);
% for i = 1:1:length(inpStruct.refData(1,:))
%     inpStruct.refData(:,i) = inpStruct.refData(:,i)-mean(inpStruct.refData(:,i));
% end

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

startInd = nlmsStruct.firOrder; % the first firOrder samples are skipped
nIter = nlmsStruct.sigLen*inpStruct.fSamp;
endInd = startInd+nIter-1;

itNo = 1;

err = zeros(nIter,1);
xEst = zeros(nIter,1);
for i = startInd:1:endInd
    smallMat = flipud(inpStruct.refData(itNo:i,:));
    xEst(itNo,1) = sum((nlmsStruct.filt).*smallMat,"all");
    err(itNo,1) = inpStruct.tarData(i,1)-xEst(itNo,1);
    ssT = sum(smallMat.^2,"all");
    ssTNorm = ssT+nlmsStruct.delta; % avoid singularity
    corrMat = nlmsStruct.alpha/ssTNorm*(smallMat*err(itNo,1));
    
    nlmsStruct.filt = nlmsStruct.filt + corrMat;
    itNo = itNo+1;
%     if(rem(itNo,10000)==0)
%         disp('Stop');
%     end
end
filtOut = nlmsStruct.filt;
disp('NLMS executed');

end