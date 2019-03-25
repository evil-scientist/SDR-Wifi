%Generate and transmit a non-HT PSDU
cfg = wlanNonHTConfig;
txPSDU = randi([0 1],8*cfg.PSDULength,1);
txSig = wlanNonHTData(txPSDU,cfg);

%Generate an L-LTF for channel estimation
txLLTF = wlanLLTF(cfg);

%Create an 802.11g channel with a 3 Hz maximum Doppler shift and a 100 ns RMS path delay. 
%Disable the reset before filtering option so that the L-LTF and data fields use the same channel realization.
ch802 = comm.RayleighChannel('SampleRate',20e6,'MaximumDopplerShift',3,'PathDelays',100e-9);

%Pass the L-LTF and data signals through an 802.11g channel with AWGN.
rxLLTF = awgn(ch802(txLLTF),10);
rxSig = awgn(ch802(txSig),10);

%Demodulate the L-LTF and use it to estimate the fading channel.
dLLTF = wlanLLTFDemodulate(rxLLTF,cfg);
chEst = JS_wlanLLTFChannelEstimateSpecial(dLLTF,cfg);

%Recover the non-HT data using the L-LTF channel estimate 
%and determine the number of bit errors in the transmitted packet.
rxPSDU = wlanNonHTDataRecover(rxSig,chEst,0.1,cfg);
[numErr,ber] = biterr(txPSDU,rxPSDU)