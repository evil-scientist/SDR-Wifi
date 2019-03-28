%% 802.11n (HT) Packet Error Rate Simulation for 2x2 TGn Channel (fading with AWGN)
%% Processing SNR Points
% For each SNR point a number of packets are tested and the packet error
% rate calculated.
%
% For each packet the following processing steps occur:
%
% 1) A PSDU is created and encoded to create a single packet waveform.

% 2) The waveform is passed through a different realization of the TGn
% channel model.

% 3) AWGN is added to the received waveform to create the desired average
% SNR per subcarrier after OFDM demodulation.
% 'comm.AWGNChannel' is configured to provide the correct SNR. 
% The configuration accounts for normalization
% within the channel by the number of receive antennas, and the noise
% energy in unused subcarriers which are removed during OFDM demodulation.

% 4) The packet is detected.

% 5) Coarse carrier frequency offset is estimated and corrected.

% 6) Fine timing synchronization is established. The L-STF, L-LTF and L-SIG
% samples are provided for fine timing to allow for packet detection at the
% start or end of the L-STF.

% 7) Fine carrier frequency offset is estimated and corrected.

% 8) The HT-LTF is extracted from the synchronized received waveform. The
% HT-LTF is OFDM demodulated and channel estimation is performed.

% 9) The HT Data field is extracted from the synchronized received waveform.
% The PSDU is recovered using the extracted field and the channel estimate.
%%
clc;
clear all;
% Create a format configuration object for a 2-by-2 HT transmission
cfgHT = wlanHTConfig;
cfgHT.ChannelBandwidth = 'CBW20'; % 20 MHz channel bandwidth
cfgHT.NumTransmitAntennas = 1;    % transmit antennas
cfgHT.NumSpaceTimeStreams = 1;    % space-time streams
cfgHT.PSDULength = 1500;          % PSDU length in bytes
cfgHT.MCS = 7;                    % spatial streams and mod & cod scheems 64-QAM rate-5/6 
cfgHT.ChannelCoding = 'LDPC';      % BCC channel coding % no space time block coding

% Channel Configuration
% In this example a TGn N-LOS channel model is used with delay profile
% Model-B. For Model-B when the distance between transmitter and receiver
% is greater than or equal to five meters, the model is NLOS. This is
% described further in 'wlanTGnChannel'.

% Create and configure the channel
tgnChannel = wlanTGnChannel;
tgnChannel.DelayProfile = 'Model-B';
tgnChannel.NumTransmitAntennas = cfgHT.NumTransmitAntennas;
tgnChannel.NumReceiveAntennas = 1;
tgnChannel.TransmitReceiveDistance = 10; % Distance in meters for NLOS
tgnChannel.LargeScaleFadingEffect = 'None';
%tgnChannel.NormalizePathGains = false; 
fs = wlanSampleRate(cfgHT); % Get the baseband sampling rate
info(tgnChannel);
[htData,htPilots] = helperSubcarrierIndices(cfgHT,'HT'); % Get the number of occupied subcarriers in HT fields and FFT length
Nst_ht = numel(htData)+numel(htPilots);
Nfft = helperFFTLength(cfgHT);      % FFT length

tgnChannel.SampleRate = fs; % Set the sampling rate of the channel
ind = wlanFieldIndices(cfgHT); % Indices for accessing each field within the time-domain packet

snr = 10:5:35;% Simulation Parameters
S = numel(snr);

packetErrorRate_LS = zeros(S,1);
packetErrorRate_MMSE = zeros(S,1);
numBitErrorsRate_LS = zeros(S,1);
numBitErrorsRate_MMSE = zeros(S,1);
BER_LS = zeros(S,1);
BER_MMSE = zeros(S,1);
maxNumPackets = 200; % Maximum number of packets at an SNR point (1/maxNumPackets of PER resolution)
use_maxNumPEs = 0; %bit dangerous to use... look below about the resolution
if use_maxNumPEs
    maxNumPEs = 10; % The maximum number of packet errors at an SNR point %max PER != 1 is (1-maxNumPEs/MaxNumPackets)
else
    maxNumPEs = maxNumPackets;
end
%parfor i = 1:S % Use 'parfor' to speed up the simulation
for i = 1:S % Use 'for' to debug the simulation 
    stream = RandStream('combRecursive','Seed',0); % Set random substream index per iteration to ensure that each
    stream.Substream = i; % iteration uses a repeatable set of random numbers
    RandStream.setGlobalStream(stream);
    
    numPacketErrors = 0;
    missed_LSTF_Packets = 0;
    packet_offeseted = 0;
    numPacketErrors_LS = 0;
    numBitErrors_LS = 0;
    numPacketErrors_MMSE = 0;
    numBitErrors_MMSE = 0;
    bitErrorRate_LS = 0;
    bitErrorRate_MMSE = 0;

    n = 1; % Index of packet transmitted
    while n<=maxNumPackets % Loop to simulate multiple packets
        
        txPSDU = randi([0 1],cfgHT.PSDULength*8,1); % PSDULength in bytes
        txLTF = wlanHTLTF(cfgHT);
        txDataSig = wlanHTData(txPSDU,cfgHT);
        
        chNoise = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',snr(i));
        
        rxLTF = chNoise(tgnChannel(txLTF));
        rxDataSig = chNoise(tgnChannel(txDataSig));
        
        nVar = 10^((-228.6 + 10*log10(290) + 10*log10(fs) + 30)/10);
        awgnChan = comm.AWGNChannel('NoiseMethod','Variance','Variance',nVar);
        
        rxLTF = awgnChan(rxLTF);
        rxDataSig = awgnChan(rxDataSig);
        
        dLTF = wlanHTLTFDemodulate(rxLTF,cfgHT);
        [chanEst_MMSE, chanEst_LS] = wlanHTLTFChannelEstimateSury(dLTF,cfgHT,nVar);
        
        reset(tgnChannel); % Reset channel for different realization
          
        rxPSDU_MMSE = wlanHTDataRecover(rxDataSig,chanEst_MMSE,nVar,cfgHT);
        [bit_err,err] = biterr(txPSDU,rxPSDU_MMSE);
        numBitErrors_MMSE = numBitErrors_MMSE + bit_err;
        bitErrorRate_MMSE = bitErrorRate_MMSE + err; 
        packetError = any(bit_err);
        numPacketErrors_MMSE = numPacketErrors_MMSE + packetError;
       
        rxPSDU_LS = wlanHTDataRecover(rxDataSig,chanEst_LS,nVar,cfgHT);
        [bit_err,err] = biterr(txPSDU,rxPSDU_LS);
         numBitErrors_LS = numBitErrors_LS + bit_err;
         bitErrorRate_LS = bitErrorRate_LS + err; 
         packetError = any(bit_err);
         numPacketErrors_LS = numPacketErrors_LS + packetError;
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        n = n+1;
    
    end
    numBitErrorsRate_LS(i) = numBitErrors_LS/(length(txPSDU) * maxNumPackets);
    numBitErrorsRate_MMSE(i) = numBitErrors_MMSE/(length(txPSDU) * maxNumPackets);
    
    BER_LS(i) = bitErrorRate_LS / maxNumPackets;
    BER_MMSE(i) = bitErrorRate_MMSE / maxNumPackets;
    % Calculate packet error rate (PER) at SNR point for LS
    if numPacketErrors_LS < maxNumPEs
        packetErrorRate_LS(i) = numPacketErrors_LS/maxNumPackets;
    else 
        packetErrorRate_LS(i) = 1;
    end
    disp(['LS: SNR ' num2str(snr(i)) ' completed after '  num2str(n-1) ' packets,' ' PER: ' num2str(packetErrorRate_LS(i))]);
      
    % Calculate packet error rate (PER) at SNR point for MMSE
    if numPacketErrors_MMSE < maxNumPEs
        packetErrorRate_MMSE(i) = numPacketErrors_MMSE/maxNumPackets;
    else 
        packetErrorRate_MMSE(i) = 1;
    end
    disp(['MMSE: SNR ' num2str(snr(i)) ' completed after ' num2str(n-1) ' packets,' ' PER: ' num2str(packetErrorRate_MMSE(i))]);
end

%% Plot Packet Error Rate vs SNR Results
figure;
semilogy(snr,packetErrorRate_LS,'-ob', 'Color', 'r');
hold on;
semilogy(snr,packetErrorRate_MMSE,'-x', 'Color', 'b');
grid on;
xlabel('SNR [dB]');
ylabel('PER');
title(join(['802.11n ', cfgHT.ChannelBandwidth(4:end), ...
    'MHz, MCS', num2str(cfgHT.MCS), ...
    ', Direct Mapping, ', num2str(cfgHT.NumTransmitAntennas), ...
    'x', num2str(tgnChannel.NumReceiveAntennas) , ...
    ', Channel ', tgnChannel.DelayProfile]));
legend('LS', 'MMSE-DFT');


figure;
semilogy(snr,BER_LS,'-ob', 'Color', 'r');
hold on;
semilogy(snr,BER_MMSE,'-x', 'Color', 'b');
grid on;
xlabel('SNR [dB]');
ylabel('BER');
title(join(['802.11n ', cfgHT.ChannelBandwidth(4:end), ...
    'MHz, MCS', num2str(cfgHT.MCS), ...
    ', Direct Mapping, ', num2str(cfgHT.NumTransmitAntennas), ...
    'x', num2str(tgnChannel.NumReceiveAntennas) , ...
    ', Channel ', tgnChannel.DelayProfile]));
legend('LS', 'MMSE-DCT');

N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
yMean = mean(y);                                    % Mean Of All Experiments At Each Value Of ‘x’
ySEM = std(y)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

