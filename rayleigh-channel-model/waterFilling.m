function P = waterFilling(SNR,Pmax)
    % initial power allocation
    initP = (Pmax + sum(1./SNR)) ./ ( length(SNR) ) - 1./SNR;

    % waterfilling algorithm
    while any( initP < 0 )
        negIndex        = initP <= 0;
        posIndex        = initP >  0;
        NkRem           = nnz(posIndex); % # of non-zero elements in posIndex 
        SNRRem          = SNR(posIndex); 
        powAllcTemp     = (Pmax + sum(1./SNRRem)) ./ (NkRem) - 1./SNRRem;
        initP(negIndex) = 0;
        initP(posIndex) = powAllcTemp;
    end
    P              = initP;
end