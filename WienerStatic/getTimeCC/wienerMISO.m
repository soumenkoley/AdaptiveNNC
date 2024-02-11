function [Rxx,Pxy] = wienerMISO(inpStruct,firOrder)

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
    
    Rxx = zeros(2*firOrder+1,nRefCha/2*(nRefCha+1));
    Pxy = zeros(2*firOrder+1,nRefCha);
    
    ccCount = 1;
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
            Rxx(:,ccCount) = rx;
            ccCount = ccCount+1;
        end
    end

    ccCount = 1;
    for s = 1:nRefCha
        x = inpStruct.refData(:,s);
        if(inpStruct.scaleInp==1)
            px = xcorr(inpStruct.tarData-mean(inpStruct.tarData),x-mean(x),firOrder);
        else
            px = xcorr(inpStruct.tarData,x,firOrder);
        end

        Pxy(:,ccCount) = px;
        ccCount = ccCount+1;
    end
    
end