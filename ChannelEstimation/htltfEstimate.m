function [est, H_LS_est] = htltfEstimate(sym,chanBW,numSTS,numESS,ind)
%htltfEstimate Channel estimate using the HT-LTF
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   EST = htltfEstimate(SYM,CHANBW,NUMSTS,NUMESS,IND) returns channel
%   estimate for each subcarrier specified by the indices IND, using
%   received symbols SYM, channel bandwidth CHANBW, number of space-time
%   streams NUMSTS and number of extension streams NUMESS.

%   Copyright 2015-2016 The MathWorks, Inc.

%#codegen

if ((numSTS+numESS)==1)
    % If one space time stream then use LS estimation directly
    ltf = wlan.internal.vhtltfSequence(chanBW,numSTS,numESS);
    %%%%%%%%
    Y = squeeze(sym(:,1,:)); %received symbols
    X = ltf(ind); %sent symbols (select only the same carrier symbols)
    N = length(ind); %Nº of carriers = len X = len Y = len(est) = len(H_LS_est)
    %%%%%%%%  LS   %%%%%%%%%
    H_LS_est = bsxfun(@rdivide,Y,X);
    H_LS_est = permute(H_LS_est,[1 3 2]);
    %%%%%%%DFT-CE%%%%%%
    F = dftmtx(N) ./sqrt(N);
    F_herm = conj(F');
    X = diag(X); %X must be a matrix now!
    X_herm = conj(X');
    I_N = eye(N);
    Lg = 16; %nº of cyclic prefixes
    %nVar = 10^((-228.6 + 10*log10(290) + 10*log10(20e6))/10); %the 20e6 should be fs = B %try 0.006 or use helper
    h = ifft(H_LS_est);
    nVar = sum(abs(h(Lg:N)')*abs(h(Lg:N)))/(N-Lg);
    R_hh = h* conj(h');
    H_MMSE_est = F * R_hh * F_herm * X_herm * inv(X * F * R_hh * F_herm * X_herm + nVar .* I_N) * Y;
    %%%%%%%%%%%%
    est = H_MMSE_est;
else               
    % MIMO channel estimation as per Perahia, Eldad, and Robert Stacey.
    % Next Generation Wireless LANs: 802.11 n and 802.11 ac. Cambridge
    % university press, 2013, page 100, Eq 4.39.
    [ltf,P,dltf,eltf] = wlan.internal.vhtltfSequence(chanBW,numSTS,numESS);

    % Verify enough symbols to estimate
    nsym = size(sym,2);
    coder.internal.errorIf(nsym<dltf+eltf, ...
        'wlan:wlanChannelEstimate:NotEnoughHTSymbols',numSTS,numESS, ...
        dltf+eltf,nsym);

    Pd = P(1:numSTS,1:dltf)'; % Extract, conjugate P matrix for HT-DLTFs
    denomD = dltf.*ltf(ind);

    Pe = P(1:numESS,1:eltf)'; % Extract, conjugate P matrix for HT-ELTFs
    denomE = eltf.*ltf(ind);
    numRx = size(sym,3);
    
    est = complex(zeros(numel(ind),numSTS+numESS,numRx));
    for i = 1:numRx
        rxsym = squeeze(sym(:,(1:dltf),i)); % Symbols on 1 receive antenna
        for j = 1:numSTS
            est(:,j,i) = rxsym*Pd(:,j)./denomD;
        end
        
        if numESS>0
            rxsym = squeeze(sym(:,dltf+(1:eltf),i)); % Symbols on 1 receive antenna
            for j = 1:numESS
                est(:,numSTS+j,i) = rxsym*Pe(:,j)./denomE;
            end
        end
    end
end
end