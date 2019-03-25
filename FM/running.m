
%y=loadFile('FMcapture1.dat'); 

RX = comm.SDRRTLReceiver('0','CenterFrequency',104.60e6,...
'EnableTunerAGC',false,'TunerGain',20,'SampleRate',2.5e6,...
'OutputDataType','single','FrequencyCorrection',0);
y = step(RX);


d = decimate(y,8,'fir'); % decimate == LPF + Down sampling to 312.5 KHz
plot_FFT_IQ(d,1,.002*2.5E6/8,2.5/8,104.600,'Spectrum of decimated signal');
figure();
[y_FM_demodulated] = FM_IQ_Demod(d); %d is the decimated signal
plot_FFT_IQ(y_FM_demodulated,1,.05*2.5E6/8,2.5/8,0,'Spectrum of baseband signal');
figure();
df = decimate(y_FM_demodulated,10,'fir'); % decimate again to bring frequency to 31.25KHz
plot_FFT_IQ(df,1,.05*2.5E6/8/10,2.5/8/10,0,'Spectrum of demodulated signal');
sound(df,2.5E6/8/10);