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
sum_rate_loss_rss = zeros(length(SNR_values), 1);

sum_rate_zf = zeros(length(SNR_values), 1);
 sum_rate_loss_zf = zeros(length(SNR_values), 1);
 disp(['Feedback budget = ' num2str(B*2) ' bits,  BS Antennas = ' num2str(M,'%02d') ]);

for t = 1:length(P_values)
P= P_values(t);
disp([ 'SNR = ' num2str(t) ' dB']);
% sum_rate_Pzf(t)= mud_scheme_comp(B,regions,M,K,P, 'P_ZF',  'Norm',MaxIter);%zeroforcing perfect csit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mud_scheme_comp start
    BFType = 0; %perfect csit
    QuantizType=0;   %Perfect CSI
    CQIType = 0;

     temp_rate = 0;
    for its = 1:MaxIter

        H = 1/sqrt(2) * (randn(K, M) + 1i * randn(K, M)); %channel
      
%              rate_temp = compute_ZF_rate(H,H,P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start perfect compute_ZF_rate 

num_users = size(H, 1);
zf_F = (H' * inv(H * H'));
for k = 1:num_users
    zf_F(:,k) = sqrt(P/num_users) * zf_F(:,k)/norm(zf_F(:,k));  
end;
gain_matrix = abs(H * zf_F).^2;
te = 1+ sum(gain_matrix,2);
temp = diag(te) * ones(num_users);
int_power = temp - gain_matrix;
SINR_matrix = gain_matrix ./ int_power;
rate = sum(log2(1 + diag(SINR_matrix)));
rate_temp=rate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end perfect  compute_ZF_rate
 
 temp_rate=temp_rate+rate_temp;
    end; 
        sum_rate = temp_rate/MaxIter;
        
sum_rate_loss=0;
 sum_rate_Pzf(t)=sum_rate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mud_scheme_comp end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sum_rate_zf(t)= mud_scheme_comp(B,regions,M,K,P, 'ZF',  'Norm',MaxIter);%zeroforcing & RVQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mud_scheme_comp zf-RVQ start

BFType = 1; %zf
QuantizType=1;%Random vector quantization
    CQIType = 0;
     temp_rate2 = 0;
    for its = 1:MaxIter
        H = 1/sqrt(2) * (randn(K, M) + 1i * randn(K, M)); %channel
            QuantH = quantiz_channels(H, regions, P, QuantizType, CQIType);
%                          rate_temp2 = compute_ZF_rate(QuantH,H, P);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%start zero forcing
     
num_users = size(H, 1);
zf_F = (QuantH' * inv(QuantH * QuantH'));%wheight matrix
for k = 1:num_users
    zf_F(:,k) = sqrt(P/num_users) * zf_F(:,k)/norm(zf_F(:,k));%w*a
end;
gain_matrix = abs(H * zf_F).^2;% 
temp = 1+ sum(gain_matrix,2);%noies+int for user
temp = diag(temp) * ones(num_users);%all interference
int_power = temp - gain_matrix;%interference power for others
SINR_matrix = gain_matrix ./ int_power;
rate = sum(log2(1 + diag(SINR_matrix)));
   rate_temp2=rate;       
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end zero forcing

      
 temp_rate2=temp_rate2+rate_temp2;
    end; 
        sum_rate2 = temp_rate2/MaxIter;
        
% sum_rate_loss2=0;
 sum_rate_zf(t)=sum_rate2 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mud_scheme_comp zf-RVQ end
sum_rate_loss_zf(t)=sum_rate_Pzf(t)-sum_rate_zf(t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RSS


BFType = 1; %zf
QuantizType=1;%Random vector quantization
    CQIType = 0;
     temp_rate3 = 0;
    for its = 1:MaxIter
        H = 1/sqrt(2) * (randn(K, M) + 1i * randn(K, M)); %channel
            QuantH = quantiz_channels(H, regions, P, QuantizType, CQIType);
%                          rate_temp2 = compute_ZF_rate(QuantH,H, P);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start zero forcing     
num_users = size(H, 1);
%EQ 19
if B<= (M-1)*log2(P)-(M-1)*log2(2*(M-1)/M)-(M-1)*log2(exp(1)-1)
r=1/(P*M/(2*(M-1))*2^(-B/(M-1))+2-exp(1));
% r=2*(M-1)/(P*M)*2^(B/(M-1));
else
r=1;
end

Pc=P*(1-r);
Pu=P*r/num_users;
Wc=zeros(M,num_users);

zf_F = (QuantH' * inv(QuantH * QuantH'));%wheight matrix:The zero-forcing (ZF) technique nulli?es the interference by the following weight matrix:
Wc(:,1)=exp(1*2*pi*rand(4,1));
Wc(:,1)= sqrt(Pc) * Wc(:,1)/norm(Wc(:,1)); 

for k = 1:num_users
    zf_F(:,k) = sqrt(Pu) * zf_F(:,k)/norm(zf_F(:,k));  

Wc(:,k)=Wc(:,1);

end
% 
gain_matrix = abs(H * zf_F).^2;
gain_matrix2=abs(H * Wc).^2;

temp = 1+ sum(gain_matrix,2);
temp = diag(temp) * ones(num_users);
int_power = temp - gain_matrix;
% int_power = temp - gain_matrix+gain_matrix2;

SINR_matrix1 =diag( gain_matrix ./ int_power);
SINR_matrix2= diag(gain_matrix2 ./ temp);
% SINR_matrixc=min(SINR_matrix2);
SINR_matrixc=sum((SINR_matrix2));

rate1 = sum(log2(1 + SINR_matrix1));
rate2=log2(1 + SINR_matrixc);
rate_temp3=rate1+rate2;
 
       temp_rate3=temp_rate3+rate_temp3;
    end; 
         sum_rate3 = temp_rate3/MaxIter;

        
sum_rate_loss2=0;
 sum_rate_zf3(t)=sum_rate3 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mud_scheme_comp zf-RVQ end
sum_rate_loss_zf(t)=sum_rate_Pzf(t)-sum_rate_zf(t);
sum_rate_loss_rss(t)=sum_rate_Pzf(t)-sum_rate_zf3(t);

end

sum_rate_loss_theor_zf=2*log2(1+P_values.*M/(2*(M-1))*2^(-B/(M-1)));
hold on
s=[0:2:40]
c=[ 2.26 2.98 3.78 4.62 5.48 6.36... 
    7.2  7.93 8.8 9.55 10.02...
     10.67 11.36 12.1 12.82 13.35...
     13.96 14.64 15.38  16 16.97]
    

figure; hold on;
plot(SNR_values, sum_rate_Pzf, 'b');
plot(SNR_values, sum_rate_zf, 'r');
plot(SNR_values, sum_rate_zf3, 'k');
% plot(SNR_values,c(:,1:11))
xlabel('SNR'); ylabel('Sum Rate (bps/hz)');
title('Sum rate performance, M= 4, B= 10');
legend('Perfect-ZFBF','RVQ-ZFBF','Approximated RSS','GA-based RSS');


% sum_rate_loss_rss_ga=sum_rate_Pzf-c(:,1:11)';

figure; hold on;
plot(SNR_values, sum_rate_loss_zf, 'r');
plot(SNR_values, sum_rate_loss_rss, '^-k');
% plot(SNR_values, sum_rate_loss_rss_ga, '^-k');

% plot(SNR_values, sum_rate_loss_theor_zf, 'b');
xlabel('SNR'); ylabel('Sum Rate Loss(bps/hz)');
title('Sum rate loss, M= 4, B= 10');
legend('ZFBF RVQ MC','\Delta(R_s)( t_e_q^s)','GA-based');



