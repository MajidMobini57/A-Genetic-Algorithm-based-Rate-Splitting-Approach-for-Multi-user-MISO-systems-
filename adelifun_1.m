function zarfiat=adelifun_1(r)
B=10;%feedback budget bits
regions = 2^B;%for code book of quantization
K=2;%K=number of users?
M=4;% BS Antennas
MaxIter=1000;
SNR_values = [0:2:40];
P_values = 10.^(SNR_values./10);
sum_rate_Pzf = zeros(length(SNR_values), 1);
sum_rate_loss_Pzf = zeros(length(SNR_values), 1);
sum_rate_loss_rss = zeros(length(SNR_values), 1);

sum_rate_zf = zeros(length(SNR_values), 1);
 sum_rate_loss_zf = zeros(length(SNR_values), 1);
 disp(['Feedback budget = ' num2str(B*2) ' bits,  BS Antennas = ' num2str(M,'%02d') ]);

P= P_values(21);
disp([ 'SNR = ' num2str(20) ' dB']);


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
%  if B>= (M-1)*log2(P)-(M-1)*log2(2*(M-1)/M)-(M-1)*log2(exp(1)-1)
% 
% 
% %  r=1/(P*M/(2*(M-1))*2^(-B/(M-1))+2-exp(1));
%  %r=2*(M-1)/(P*M)*2^(B/(M-1));
% % else
% r=1;
%  end
pc=P*(1-r(1));
Pu=P*r(1)/num_users;
Wc=zeros(M,num_users);
zf_F = (QuantH' * inv(QuantH * QuantH'));%wheight matrix:The zero-forcing (ZF) technique nulli?es the interference by the following weight matrix:
Wc(:,1)=exp(1*2*pi*rand(4,1));
Wc(:,1)= sqrt(pc) * Wc(:,1)/norm(Wc(:,1)); 

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
 SINR_matrixc=min(SINR_matrix2);
% SINR_matrixc=sum((SINR_matrix2));

rate1 = sum(log2(1 + SINR_matrix1));
rate2=log2(1 + SINR_matrixc);
rate_temp3=rate1+rate2;
 
       temp_rate3=temp_rate3+rate_temp3;
    end; 
         sum_rate3 = temp_rate3/MaxIter;

        
sum_rate_loss2=0;
 capacity=sum_rate3 
 zarfiat=1/capacity;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mud_scheme_comp zf-RVQ end

sum_rate_loss_theor_zf=2*log2(1+P_values.*M/(2*(M-1))*2^(-B/(M-1)));
