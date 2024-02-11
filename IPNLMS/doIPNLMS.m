function [ipnlms,err, dEst, mu] = doIPNLMS(ipnlms,inpStruct,P,L)
% this function implements the IPNLMS algorithm

%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

startInd = L; % the first firOrder samples are skipped
nIter = ipnlms.sigLen*ipnlms.fSamp;
endInd = startInd+nIter-1;

dEst = zeros(nIter,1);
err = zeros(nIter,1);
mu = zeros(nIter,1);

itNo = 1;
if(ipnlms.start==1)
    ipnlms.X(1:(L-1),:) = flipud(inpStruct.refData(1:(startInd-1),:));
    ipnlms.start = 0;
end
for i = startInd:1:endInd
    ipnlms.X(2:L,:) = ipnlms.X(1:(L-1),:);
    ipnlms.X(1,1:P) = (inpStruct.refData(i,1:P)); % new data entry

    ipnlms.sigmaX = sum(ipnlms.X.^2,'all')/L;
    ipnlms.deltaNLMS = ipnlms.cst*ipnlms.sigmaX;
    %ipnlms.deltaNLMS = 10^-8;
    %if(itNo>sigmaBuffSamp)
        %startBuff = 0;
        % computation starts
    dEst(itNo,1) = sum((ipnlms.X).*(ipnlms.H),'all');
    err(itNo,1) = inpStruct.tarData(i,1) - dEst(itNo,1);
    hNorm = sum(abs(ipnlms.H),1); % 1xP vector
    
    % g is a LXP matrix
    g = (1-ipnlms.K)/(2*L) + (1+ipnlms.K)*(abs(ipnlms.H)./(2*repmat(hNorm,L,1) + ipnlms.ep));
        
    %den = sum(((ipnlms.X).^2).*g) + ipnlms.cst*sigmaAvgUse*(1-ipnlms.K)/(2*PL);
    sumX = sum(((ipnlms.X).^2).*g,'all');
    
    den = sumX + (1-ipnlms.K)/(2*L)*ipnlms.deltaNLMS;
    
    mu = ipnlms.alpha/den; % LXP matrix
    ipnlms.H = ipnlms.H + mu.*(err(itNo,1)*(g.*ipnlms.X));
    %end

%     if(rem(itNo,10000)==0)
%         disp('Stop');
%     end
    itNo = itNo + 1;
end

end
