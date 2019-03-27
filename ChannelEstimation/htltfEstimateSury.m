function [est, H_LS_est] = htltfEstimate2(sym,chanBW,numSTS,numESS,ind, snr)
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
    %F_herm = conj(F');
    F_herm = F';
    X = diag(X); %X must be a matrix now!
    %X_herm = conj(X');
    X_herm = X';
    I_N = eye(N);
    Lg = 16; %nº of cyclic prefixes %or 16?
    h_pre = ifft(H_LS_est);
    %disp(h);
    nVar = 0;
    for foo = Lg+1:N
        nVar = nVar + abs(h_pre(foo));
        %abs(h(foo));
    end
    nVar = nVar/(N-Lg);
    %disp('foo');
    %disp(nVar);
    %disp('Notfoo');
    %disp(sum(abs(h(Lg+1:N)')*abs(h(Lg+1:N)))/(N-Lg));
    h = zeros(1,56);
    for foo = 1:Lg
        h(foo) = h_pre(foo);
    end
    %R_hh = h*h';
    %disp(length(h));
    [xc,lags] = xcorr(h,h,length(h)-1);
    r = xc(length(h):end);
    %r = xc(1:length(h));
    %disp(size(r));
    R_hh = toeplitz(r,conj(r));
    %disp(F * R_hh * F_herm * X_herm);
    %disp(size(R_hh));    
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
    
    H_LS_est = complex(zeros(numel(ind),numSTS+numESS,numRx));
    for i = 1:numRx
        rxsym = squeeze(sym(:,(1:dltf),i)); % Symbols on 1 receive antenna
        for j = 1:numSTS
            H_LS_est(:,j,i) = rxsym*Pd(:,j)./denomD;
        end
    end
    %disp(size(H_LS_est));
    %%%%%%%%%%%%%%%%   
    %disp(ltf(ind)');
    H_tilda = complex(zeros(numel(ltf),numSTS+numESS,numRx));
    h_tilda = complex(zeros(numel(ltf),numSTS+numESS,numRx));
    %disp(size(H_tilda));
    for i = 1:numRx
        for j = 1:numSTS
            %for k = 1:numel(ltf)
                %if k == 0 
                    H_tilda(1,j,i) = 0;
                %elseif k < 29
                    H_tilda(2:29,j,i) = H_LS_est(1:28,j,i);
                %elseif k < 36
                    H_tilda(30:36,j,i) = 0;
                %elseif k < 64
                    H_tilda(37:64,j,i) = H_LS_est(29:56,j,i);
                %end
            %end
        end
    end
    %disp(size(H_tilda));
    for i = 1:numRx
        for j = 1:numSTS
            h_tilda(:,j,i) = ifft(H_LS_est(:,j,i));
        end
    end
    
    for i = 1:numRx
        for j = 1:numSTS
            for k = 1:numel(ltf)
                if k > 16 % RMS Delay Spread of Channel for B == 15
                    h_tilda(k,j,i) = 0;
                end
            end
        end
    end
%    disp(h_tilda);
    
    est_64 = complex(zeros(numel(ltf),numSTS+numESS,numRx));
    for i = 1:numRx
        for j = 1:numSTS
            est_64(:,j,i) = fft(h_tilda(:,j,i));
        end    
    end
    
    %est_26 = est_64(5:30,:,:);
    %est_27 = est_64(32:61,:,:);
    %est = [est_26 ; est_27];
    est = est_64(1:56,:,:);
    
    %disp(size(est));
 end
end