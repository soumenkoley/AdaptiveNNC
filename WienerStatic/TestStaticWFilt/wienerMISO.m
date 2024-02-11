function [hh,rankOut] = wienerMISO(inpStruct,firOrder)

    % this function is a Wiener filter implementation for the 
    % MISO system
    % https://pubs.geoscienceworld.org/geophysics/article-abstract/35/5/785/71204/principles-of-digital-multichannel-filtering

    % code very similar to:
    % https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwiBpPrd-fr8AhXu_7sIHWYgCx4QFnoECAkQAQ&url=https%3A%2F%2Flabcit.ligo.caltech.edu%2F~cit40m%2FDocs%2FT070192_KeenanPepper.pdf&usg=AOvVaw2eVnTxFKaUQ9SZVvqwh5mT
    % K Pepper, 17-08-2017, LIGO-T070192-00-R
    
    %----------------------------------------------------------------------
    % author: S. Koley
    % Department of Physics
    % Gran Sasso Science Institute
    % soumen.koley@gssi.it
    % ---------------------------------------------------------------------

    [~,nRefCha] = size(inpStruct.refData);
    % problem to be formulated as Rh = P
    % R is the reference channel cross-correlation matrix
    % P is the reference to target channel cross-correlation
    % h is the set of FIR Wiener coefficients
    
    %----------------------------------------------------------------------
    % author: S. Koley
    % Department of Physics
    % Gran Sasso Science Institute
    % soumen.koley@gssi.it
    % ---------------------------------------------------------------------
    
    R = zeros(nRefCha*(firOrder+1), nRefCha*(firOrder+1));

    for s1 = 1:1:nRefCha
        disp(['Estimating correlations for reference sensor ',num2str(s1)]);
        for s2 = s1:nRefCha
            x = inpStruct.refData(:,s1);
            y = inpStruct.refData(:,s2);
            if(inpStruct.scaleInp==1)
                rx = xcorr(x-mean(x), y-mean(y), firOrder);
            else
                rx = xcorr(x, y, firOrder);
            end
            Rx = toeplitz(flipud(rx(1:(firOrder+1),1)), rx((firOrder+1):(2*firOrder+1),1));

            top = (s1-1)*(firOrder+1)+1;
            bottom = s1* (firOrder+1);
            left = (s2-1)*(firOrder+1)+1;
            right = s2*(firOrder + 1);

            R(top:bottom, left:right) = Rx;

            if(s1~=s2)
                R(left:right, top:bottom) = Rx'; % transpose
            end
        end
    end

    P = zeros(nRefCha*(firOrder+1),1);
    for s = 1:nRefCha
        x = inpStruct.refData(:,s);
        if(inpStruct.scaleInp==1)
            px = xcorr(inpStruct.tarData-mean(inpStruct.tarData),x-mean(x),firOrder);
        else
            px = xcorr(inpStruct.tarData,x,firOrder);
        end

        top = (s-1)*(firOrder + 1) + 1;
        bottom = s*(firOrder+1);

        P(top:bottom,1) = px((firOrder+1):(2*firOrder+1),1);
    end
    rankOut = rank(R');
    h = lsqminnorm(R',P);
    % QR decomposition will yield teh sam result
    %[Q,MM] = qr(R',0); QR decomposition will also give the same result
    %gg = Q*(MM'\P);
    hh = reshape(h,firOrder+1,nRefCha);
    %disp('done');
    
end