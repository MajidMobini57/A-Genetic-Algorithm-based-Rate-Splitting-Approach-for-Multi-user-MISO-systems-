clc;
close all;
clear all;
B=10;%feedback budget bits
regions = 2^B;%for code book of quantization
K=2;%K=number of users?
M=4;% BS Antennas
MaxIter=10000;
SNR_values = [0:2:40];
P_values = 10.^(SNR_values./10);
sum_rate_Pzf = zeros(length(SNR_values), 1);
sum_rate_loss_Pzf = zeros(length(SNR_values), 1);
sum_rate_zf = zeros(length(SNR_values), 1);
 sum_rate_loss_zf = zeros(length(SNR_values), 1);
 sum_rate_P1  = zeros(length(SNR_values), 1);
 sum_rate_loss_P1 = zeros(length(SNR_values), 1);
 sum_rate_loss_theor_P1 = zeros(length(SNR_values), 1);
rate_P1 = zeros(length(SNR_values), 1);
disp(['Feedback budget = ' num2str(B) ' bits,  BS Antennas = ' num2str(M,'%02d') ]);

for t = 1:length(P_values)
    P= P_values(t);
disp([ 'SNR = ' num2str(t*5) ' dB']);
sum_rate_Pzf(t)= mud_scheme_comp(B,regions,M,K,P, 'P_ZF',  'Norm',MaxIter);%zeroforcing perfect csit
sum_rate_zf(t)= mud_scheme_comp(B,regions,M,K,P, 'ZF',  'Norm',MaxIter);%zeroforcing & RVQ
sum_rate_P1(t)= mud_scheme_comp(B,regions,M,K,P, 'P1',  'Norm',MaxIter);%Single-user beamforming & RVQ

sum_rate_loss_zf(t)=sum_rate_Pzf(t)-sum_rate_zf(t);
sum_rate_loss_P1(t)=sum_rate_Pzf(t)-sum_rate_P1(t);

if B<= (M-1)*log2(P)-(M-1)*log2(2*(M-1)/M)-(M-1)*log2(exp(1)-1)
sum_rate_loss_theor_P1(t)=log2(exp(1))+log2(P.*M/(M-1)*2^(-B/(M-1))+2-exp(1));
else
sum_rate_loss_theor_P1(t)=2*log2(1+P.*M/(2*(M-1))*2^(-B/(M-1)));
end
 
end
sum_rate_loss_theor_zf=2*log2(1+P_values.*M/(2*(M-1))*2^(-B/(M-1)));


figure; hold on;
plot(SNR_values, sum_rate_Pzf, 'b');
plot(SNR_values, sum_rate_zf, 'r');
% plot(SNR_values, sum_rate_P1, 'Y');

xlabel('SNR'); ylabel('Sum Rate (bps/hz)');
title('Sum rate performance, M= 4, B= 10');
% legend('Perfect-ZFBF','RVQ-ZFBF', 'P1');
legend('Perfect-ZFBF','RVQ-ZFBF');

figure; hold on;

plot(SNR_values, sum_rate_loss_zf, 'r');
% plot(SNR_values, sum_rate_loss_P1, 'Y');
plot(SNR_values, sum_rate_loss_theor_zf, 'b');
% plot(SNR_values, sum_rate_loss_theor_P1, 'g');
xlabel('SNR'); ylabel('Sum Rate Loss(bps/hz)');
title('Sum rate loss, M= 4, B= 10');
% legend('ZFBF Monte Carlo','P1 Monte Carlo', 'ZFBF Theoretical','P1 Theoretical');
legend('ZFBF Monte Carlo',                    'ZFBF Theoretical'                 );
%



