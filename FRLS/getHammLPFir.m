function [h] = getHammLPFir(M,fCut,fSamp)
    %this function was written to derive the hamming window
    %filter coefficients
    %for detailed derivation refer to my notebook or
    %watch https://www.youtube.com/watch?v=dQPEiDg1ZJo&ab_channel=EnggClasses
    
    % M = filter order
    % fCut = low-pass frequency in Hz
    % fSamp = sampling frequency in Hz

    %----------------------------------------------------------------------
    % author: S. Koley
    % Department of Physics
    % Gran Sasso Science Institute
    % soumen.koley@gssi.it
    % ---------------------------------------------------------------------

    if rem(M,2) ~= 0
        tau = (M-1)/2;
    else
        error('Filter must have an odd order');
    end

    %now get the fCut as fraction of pi where pi = Fs/2
    K = fSamp/(2*fCut);

    h = zeros(M,1);
    n = 0;
    for i = 0:(floor((M-1)/2)-1) % we can do with half since it is symmetric
        hammWin = 0.54 - 0.46*cos(2*pi*n/(M-1));
        h(i+1,1) = sin(pi/K*(n-tau))/(pi*(n-tau))*hammWin;
        n = n+1;
    end

    h(i+2,1) = 1/K;

    %use symmetry
    n = 1;
    for i = (floor((M-1)/2)+1):(M-1)
        h(i+1,1) = h(floor((M-1)/2)+1-n,1);
        n = n+1;
    end
end