function est = JS_wlanLLTFChannelEstimate(rxSym,cfgFormat)

% Channel bandwidth parameterized using object
cbw = cfgFormat.ChannelBandwidth;
numSC = size(rxSym,1);
numRxAnts = size(rxSym,3);
num20 = 1;

lltf = lltfReference(num20); % Get reference subcarriers

% Verify number of subcarriers to estimate
coder.internal.errorIf(numSC~=numel(lltf),'wlan:wlanChannelEstimate:IncorrectNumSC',numel(lltf),numSC);

ls = rxSym./repmat(lltf,1,size(rxSym,2),numRxAnts); % Least-square estimate   
est = mean(ls,2); % Average over the symbols
end

function ref = lltfReference(num20MHz)
    % 20 MHz reference
    [lltfLower, lltfUpper] = JS_LLTF_Sequence();

    % Replicate over number of 20 MHz segments ignoring the DC and reshape
    ref = reshape([lltfLower; lltfUpper]*ones(1,num20MHz),[],1);
end
