clc;
x = snr; 

yB = BER_MMSE_B; % Create Dependent Variable ‘Experiments’ Data
N = size(yB,1); % Number of ‘Experiments == PACKETS’ In Data Set
yMeanB = mean(yB); % Mean Of All Experiments At Each Value Of ‘x’
ySEMB = std(yB)/sqrt(N); % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI95B = tinv([0.025 0.975], N-1); % Calculate 95% Probability Intervals Of t-Distribution
yCI95B = bsxfun(@times, ySEMB, CI95B(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

yD = BER_MMSE_D; % Create Dependent Variable ‘Experiments’ Data
N = size(yD,1); % Number of ‘Experiments == PACKETS’ In Data Set
yMeanD = mean(yD); % Mean Of All Experiments At Each Value Of ‘x’
ySEMD = std(yD)/sqrt(N); % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI95D = tinv([0.025 0.975], N-1); % Calculate 95% Probability Intervals Of t-Distribution
yCI95D = bsxfun(@times, ySEMD, CI95D(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

yE = BER_MMSE_E; % Create Dependent Variable ‘Experiments’ Data
N = size(yE,1); % Number of ‘Experiments == PACKETS’ In Data Set
yMeanE = mean(yE); % Mean Of All Experiments At Each Value Of ‘x’
ySEME = std(yE)/sqrt(N); % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI95E = tinv([0.025 0.975], N-1); % Calculate 95% Probability Intervals Of t-Distribution
yCI95E = bsxfun(@times, ySEME, CI95E(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

figure;

h1 = errorbar(x,yMeanB,yCI95B(2,:),'-s','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
set(get(h1,'Parent'), 'YScale', 'log')
hold on

h2 = errorbar(x,yMeanD,yCI95D(2,:),'-s','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
set(get(h2,'Parent'), 'YScale', 'log')
hold on

h3 = errorbar(x,yMeanE,yCI95E(2,:),'-s','MarkerSize',5,'MarkerEdgeColor','yellow','MarkerFaceColor','yellow');
set(get(h3,'Parent'), 'YScale', 'log')
hold on

grid on;
xlabel('SNR [dB]');
ylabel('BER');
title(join(['802.11n ', '20 MHz, MCS 7']));

legend('Model B','Model D','Model E');