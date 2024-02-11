function [frls,dEst,epsp,gammaStore,rescue] = doSFTF(frls,inpStruct,P,L)
% this script implements the SFTF algorithm
% for the MISO system
%----------------------------------------------------------------------
% author: S. Koley
% Department of Physics
% Gran Sasso Science Institute
% soumen.koley@gssi.it
% ---------------------------------------------------------------------

PL = P*L;
dEst = zeros(frls.nSamp,1);
epsp = zeros(frls.nSamp,1);
gammaStore = zeros(frls.nSamp,1);

startInd = 101;
endInd = startInd + frls.nSamp-1;

itNo = 1;
for i = startInd:1:endInd
    
    % computation starts here
    frls.X((P+1):(PL+P),1) = frls.X(1:PL,1);
    frls.X(1:P,1) = (inpStruct.refData(i,1:P))';

    frls.eP(1:P,1) = frls.A(:,(P+1):(PL+P))*frls.X((P+1):(PL+P));
    frls.eP(1:P,1) = frls.eP(1:P,1) + frls.X(1:P,1);

    frls.CN1(1:P,1) = -frls.lambda1*frls.alpha1*frls.eP;

    frls.CN1((P+1):(PL+P),1) = frls.C + ((frls.CN1(1:P,1))'*frls.A(:,(P+1):(PL+P)))';
    
    % update
    gamma1N1 = frls.gamma1 - (frls.CN1(1:P,1))'*frls.eP;

    eN(1:P,1) = frls.eP*frls.gamma;

    frls.A(:,(P+1):(PL+P)) = frls.A(:,(P+1):(PL+P)) + eN*frls.C';

    frls.alpha1 = frls.lambda1*frls.alpha1 - (frls.CN1(1:P,1)*(frls.CN1(1:P,1))')/gamma1N1;

    CN1Ns(1:P,1) = frls.CN1((PL+1):(PL+P),1);
    
    rps(1:P,1) = -frls.lambda*frls.beta*CN1Ns;

    rpf(1:P,1) = frls.B*frls.X;

    y1(1:P,1) = rpf - rps;

    rp1(1:P,1) = rps + frls.K(1)*y1;

    rp2(1:P,1) = rps + frls.K(2)*y1;

    rp5(1:P,1) = rps + frls.K(5)*y1;

    CN1Nf(1:P,1) = -frls.lambda1*(frls.betaInv*rpf);

    CN1N(1:P,1) = CN1Ns + frls.K(4)*(CN1Nf-CN1Ns);

    frls.C(1:PL,1) = frls.CN1(1:PL,1) - (CN1N'*frls.B(:,1:(PL)))';

    gamma1s = gamma1N1 + CN1Ns'*rp5;

    CN1N2(1:P,1) = frls.K(2)*CN1Nf + (1-frls.K(2))*CN1Ns;

    gamma2N1 = gamma1s - CN1N2'*rp2;

    frls.betaInv = frls.lambda1*frls.betaInv - (CN1N2*CN1N2')/gamma2N1;

    gammas = 1/gamma1s;
    
    r1 = rp1*gammas;

    r2 = rp2*gammas;

    frls.B(:,1:(PL)) = frls.B(:,1:(PL)) + r1*frls.C';
    
    frls.beta = frls.lambda*frls.beta + r2*rp2';

    frls.gamma = frls.lambdaN*det(frls.beta)*det(frls.alpha1);

    frls.gamma1 = 1/frls.gamma;

    frls.XX((P+1):(PL+P),1) = frls.XX(1:(PL),1);
    frls.XX(1:P,1) = (inpStruct.refData(i,1:P))';
    
    gammaStore(itNo,1)= frls.gamma;

    epsp(itNo,1) = inpStruct.tarData(i,1);

    dEst(itNo,1) = frls.W'*frls.XX(1:PL,1);

    epsp(itNo,1) = epsp(itNo,1) + dEst(itNo,1);

    eps = epsp(itNo,1)*frls.gamma;

    frls.W = frls.W + eps*frls.C;
    
    % soft-constrained rescue
    if((frls.gamma>=1.0)||(frls.gamma<=0))
        rescue = 1;
        if(frls.doRescue)
            frls.mu = sum(inpStruct.refData((i-L+1):i,:).^2,'all');
            frls.alpha1 = diag(1/(frls.lambdaN*frls.mu)).*eye(P);
            
            frls.beta = (inv(frls.alpha1));
            
            frls.betaInv = inv(frls.beta);
            frls.gamma1 = 1.0;
            frls.gamma = 1.0;
            frls.X(1:(PL),1) = 0;
            frls.A(:,(P+1):(PL+P)) = 0;
            frls.B(:,1:PL) = 0;
            frls.C(1:(PL)) = 0;
            
        end
    else
        rescue = 0;
    end

    % some arbitrary checkpoint whle debugging
%     if(itNo==76421)
%         disp('Stop');
%     end


    itNo = itNo+1;

end

end