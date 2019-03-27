function [bits, mpduLength, fc] = create_beacon()
% https://nl.mathworks.com/help/wlan/ref/wlanmacframe.html

%%
SSID = 'TEST_BEACON';   % Network SSID
beaconInterval = 1;   % In Time Units (TU) => *1024 microsec (default is 100, but we are in a hurry)
TU = 1024e-6;           % Defined in standards
band = 2.4;               % 5 Ghz or 2.4 Ghz
chNum = 7;             % Channel number, 52 corresponds to 5260 MHz
                        % check the wireshark book for confirmation
bitsPerByte = 8;        % ...just to be more readable more at the end

%1º- create configuration for the beacon frame-body (what goes on the body)
%2º- create configuration for the beacon (structure of the frame,ect...)
%3º- with the configs, create the beacon frame (will give a bunch of Bytes)
%4º- convert these Bytes to bits (preparation for the modulation)

%% 1º Create Beacon frame-body configuration object
frameBodyConfig = wlanMACManagementConfig;          % Creates 1 management frame configuration object
frameBodyConfig.BeaconInterval = beaconInterval;    % Set Beacon Interval
frameBodyConfig.SSID = SSID;                        % Set name of the Network
dsElementID = 3;                                    % DS Parameter IE element ID
dsInformation = dec2hex(chNum, 2);                  % DS Parameter IE information (should be in hex)
frameBodyConfig = frameBodyConfig.addIE(dsElementID, dsInformation);  % Add DS Parameter IE to the management frame configuration obj


%% 2º Create Beacon frame configuration object
beaconFrameConfig = wlanMACFrameConfig('FrameType', 'Beacon');  %this creates the default beacon
beaconFrameConfig.ManagementConfig = frameBodyConfig;           %Append the management frame configuration obj

%% 3º Generate Beacon frame
[beacon, mpduLength] = wlanMACFrame(beaconFrameConfig);

%% 4º Convert the mpdu bytes in hexa-decimal format to bits
beacon_dec = hex2dec(beacon);
bits = reshape(de2bi(beacon_dec, bitsPerByte)', [], 1);


% Calculate center frequency for the given band and channel number
fc = helperWLANChannelFrequency(chNum, band);