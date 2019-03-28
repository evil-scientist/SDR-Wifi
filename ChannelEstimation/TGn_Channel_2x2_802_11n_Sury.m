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
tgnChannel.DelayProfile = 'Model-E';
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

snr = 10:2:35;% Simulation Parameters
S = numel(snr);

packetErrorRate_MMSE = zeros(S,1);
numBitErrorsRate_MMSE = zeros(S,1);
maxNumPackets = 1000; % Maximum number of packets at an SNR point (1/maxNumPackets of PER resolution)
BER_MMSE = zeros(maxNumPackets,S);
use_maxNumPEs = 0; %bit dangerous to use... look below about the resolution
if use_maxNumPEs
    maxNumPEs = 10; % The maximum number of packet errors at an SNR point %max PER != 1 is (1-maxNumPEs/MaxNumPackets)
else
    maxNumPEs = maxNumPackets;
end
DEBUG = 0;
%parfor i = 1:S % Use 'parfor' to speed up the simulation
for i = 1:S % Use 'for' to debug the simulation 
    stream = RandStream('combRecursive','Seed',0); % Set random substream index per iteration to ensure that each
    stream.Substream = i; % iteration uses a repeatable set of random numbers
    RandStream.setGlobalStream(stream);
    
    awgnChannel = comm.AWGNChannel; % Create an instance of the AWGN channel per SNR point simulated
    awgnChannel.NoiseMethod = 'Signal to noise ratio (SNR)';
    awgnChannel.SignalPower = 1/tgnChannel.NumReceiveAntennas;% Normalization
    awgnChannel.SNR = snr(i)-10*log10(Nfft/Nst_ht); % Account for energy in nulls
   
    numPacketErrors = 0;
    missed_LSTF_Packets = 0;
    packet_offeseted = 0;
    numPacketErrors_MMSE = 0;
    numBitErrors_MMSE = 0;
    bitErrorRate_MMSE = 0;

    n = 1; % Index of packet transmitted
    while n<=maxNumPackets % Loop to simulate multiple packets
        
        txPSDU = randi([0 1],cfgHT.PSDULength*8,1); % PSDULength in bytes
        tx = wlanWaveformGenerator(txPSDU,cfgHT); % Generate a packet waveform
        tx = [tx; zeros(15,cfgHT.NumTransmitAntennas)]; % Add trailing zeros to allow for channel filter dela
        %debug = info(tgnChannel);
        reset(tgnChannel); % Reset channel for different realization
        rx_1 = tgnChannel(tx); % Pass the waveform through the TGn channel model 
        rx = awgnChannel(rx_1); % Add noise
       
        % Packet detect and determine coarse packet offset
        coarsePktOffset = wlanPacketDetect(rx,cfgHT.ChannelBandwidth);
        if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
            numPacketErrors_MMSE = numPacketErrors_MMSE+1;
            missed_LSTF_Packets = missed_LSTF_Packets +1;
            n = n+1; 
            continue; % Go to next loop iteration
        end
        
        % Extract L-STF and perform coarse frequency offset correction
        lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:); 
        coarseFreqOff = wlanCoarseCFOEstimate(lstf,cfgHT.ChannelBandwidth);
        rx = helperFrequencyOffset(rx,fs,-coarseFreqOff);
        
        % Extract the non-HT fields and determine fine packet offset
        nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:); 
        finePktOffset = wlanSymbolTimingEstimate(nonhtfields,cfgHT.ChannelBandwidth);
        
        % Determine final packet offset
        pktOffset = coarsePktOffset+finePktOffset;
        
        % If packet detected outwith the range of expected delays from the
        % channel modeling; packet error
        if pktOffset>7
            numPacketErrors_MMSE = numPacketErrors_MMSE+1;
            packet_offeseted = packet_offeseted + 1;
            n = n+1;
            continue; % Go to next loop iteration
        end

        % Extract L-LTF and perform fine frequency offset correction
        lltf = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:); 
        fineFreqOff = wlanFineCFOEstimate(lltf,cfgHT.ChannelBandwidth);
        rx = helperFrequencyOffset(rx,fs,-fineFreqOff);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extract HT-LTF samples from the waveform, demodulate and perform
        % channel estimation
        htltf = rx(pktOffset+(ind.HTLTF(1):ind.HTLTF(2)),:);
        
        % Extract HT Data samples from the waveform
        htdata = rx(pktOffset+(ind.HTData(1):ind.HTData(2)),:);
        
        htltfDemod = wlanHTLTFDemodulate(htltf,cfgHT);
        [chanEst_MMSE, chanEst_LS] = wlanHTLTFChannelEstimateSury(htltfDemod,cfgHT,htdata);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MMSE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if numPacketErrors_MMSE < maxNumPEs
            % Estimate the noise power in HT data field for MMSE
            nVarHT_MMSE = htNoiseEstimate(htdata,chanEst_MMSE,cfgHT);

            % Recover the transmitted PSDU in HT Data for MMSE
            rxPSDU_MMSE = wlanHTDataRecover(htdata,chanEst_MMSE,nVarHT_MMSE,cfgHT);

            % Determine if any bits are in error, i.e. a packet error for
            % MMSE
            [bit_err,err] = biterr(txPSDU,rxPSDU_MMSE);
            %numBitErrors_MMSE = numBitErrors_MMSE + bit_err;
            bitErrorRate_MMSE = bitErrorRate_MMSE + err; 
            BER_MMSE(n,i) = err;
            if bit_err > 0 && DEBUG
                disp(join(['MMSE, packet nº ', num2str(n), ' has ',num2str(bit_err), ...
                    ' wrong bits.Total errors so far: ', num2str(numBitErrors_MMSE)]));
            end
            packetError = any(bit_err);
            numPacketErrors_MMSE = numPacketErrors_MMSE + packetError;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        n = n+1;
    end
    %numBitErrorsRate_MMSE(i) = numBitErrors_MMSE/(length(txPSDU) * maxNumPackets);
    %BER_MMSE(i) = bitErrorRate_MMSE / maxNumPackets;
    % Calculate packet error rate (PER) at SNR point for MMSE
    if numPacketErrors_MMSE < maxNumPEs
        packetErrorRate_MMSE(i) = numPacketErrors_MMSE/maxNumPackets;
    else 
        packetErrorRate_MMSE(i) = 1;
    end
    if DEBUG
        disp(['MMSE: SNR ' num2str(snr(i)) ' completed after ' num2str(n-1) ' packets,' ' PER: ' num2str(packetErrorRate_MMSE(i))]);
    end
end
%% Plot Packet Error Rate vs SNR Results
figure;
hold on;
semilogy(snr,packetErrorRate_MMSE,'-x', 'Color', 'b');
grid on;
xlabel('SNR [dB]');
ylabel('PER');
title(join(['802.11n ', cfgHT.ChannelBandwidth(4:end), ...
    'MHz, MCS', num2str(cfgHT.MCS), ...
    ', Direct Mapping, ', num2str(cfgHT.NumTransmitAntennas), ...
    'x', num2str(tgnChannel.NumReceiveAntennas)]));
%legend(tgnChannel.DelayProfile);

%%
%legend('LS', 'MMSE-DCT');
x = snr; 
y = BER_MMSE; % Create Dependent Variable ‘Experiments’ Data
N = size(y,1); % Number of ‘Experiments == PACKETS’ In Data Set
yMean = mean(y); % Mean Of All Experiments At Each Value Of ‘x’
ySEM = std(y)/sqrt(N); % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI95 = tinv([0.025 0.975], N-1); % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
figure;
semilogy(x,yMean,'-x', 'Color', 'b');% Plot Mean Of All Experiments
xlabel('SNR [dB]');
ylabel('BER');
hold on;
semilogy(x,yCI95+yMean,'-o', 'Color', 'r'); % Plot 95% Confidence Intervals Of All Experiments
hold off;
grid on;
title(join(['802.11n ', cfgHT.ChannelBandwidth(4:end), ...
    'MHz, MCS', num2str(cfgHT.MCS), ...
    ', Direct Mapping, ', num2str(cfgHT.NumTransmitAntennas), ...
    'x', num2str(tgnChannel.NumReceiveAntennas)]));
legend('MMSE-DCT');
