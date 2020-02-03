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
 disp(['Feedback budget = ' num2str(B) ' bits,  BS Antennas = ' num2str(M,'%02d') ]);

for t = 1:length(P_values)
P= P_values(t);
disp([ 'SNR = ' num2str(t*2) ' dB']);
% sum_rate_Pzf(t)= mud_scheme_comp(B,regions,M,K,P, 'P_ZF',  'Norm',MaxIter);%zeroforcing perfect csit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mud_scheme_comp start
    BFType = 0; %perfect csit
    QuantizType=0;   %Perfect CSI


    CQIType = 0;



     temp_rate = 0;
    for its = 1:MaxIter

        H = 1/sqrt(2) * (randn(K, M) + 1i * randn(K, M)); %channel
      
             rate_temp = compute_ZF_rate(H,H,P);


 temp_rate=temp_rate+rate_temp;
    end; 
        sum_rate = temp_rate/MaxIter;
        
sum_rate_loss=0;
 sum_rate_Pzf(t)=sum_rate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mud_scheme_comp end

% sum_rate_zf(t)= mud_scheme_comp(B,regions,M,K,P, 'ZF',  'Norm',MaxIter);%zeroforcing & RVQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mud_scheme_comp zf-RVQ start

BFType = 1; %RBF/PU2RC
QuantizType=1;%Random vector quantization
    CQIType = 0;

     temp_rate2 = 0;
    for its = 1:MaxIter
        H = 1/sqrt(2) * (randn(K, M) + 1i * randn(K, M)); %channel
            QuantH = quantiz_channels(H, regions, P, QuantizType, CQIType);
             rate_temp2 = compute_ZF_rate(QuantH,H, P);
      
 temp_rate2=temp_rate2+rate_temp2;
    end; 
        sum_rate2 = temp_rate2/MaxIter;
        
sum_rate_loss2=0;
 sum_rate_zf(t)=sum_rate2 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mud_scheme_comp zf-RVQ end

sum_rate_loss_zf(t)=sum_rate_Pzf(t)-sum_rate_zf(t);


 
end

sum_rate_loss_theor_zf=2*log2(1+P_values.*M/(2*(M-1))*2^(-B/(M-1)));
figure; hold on;
plot(SNR_values, sum_rate_Pzf, 'b');
plot(SNR_values, sum_rate_zf, 'r');
xlabel('SNR'); ylabel('Sum Rate (bps/hz)');
title('Sum rate performance, M= 4, B= 10');
legend('Perfect-ZFBF','RVQ-ZFBF');


figure; hold on;
plot(SNR_values, sum_rate_loss_zf, 'r');
plot(SNR_values, sum_rate_loss_theor_zf, 'b');
xlabel('SNR'); ylabel('Sum Rate Loss(bps/hz)');
title('Sum rate loss, M= 4, B= 10');
legend('ZFBF Monte Carlo',                    'ZFBF Theoretical'                 );
%



