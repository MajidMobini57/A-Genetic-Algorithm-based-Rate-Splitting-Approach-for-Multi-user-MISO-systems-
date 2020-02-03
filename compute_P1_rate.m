function [rate] = compute_P1_rate(QH,H,P,B,M)


num_users = size(H, 1);

if B<= (M-1)*log2(P)-(M-1)*log2(2*(M-1)/M)-(M-1)*log2(exp(1)-1)
% r=1/(P*M/(2*(M-1))*2^(-B/(M-1))+2-exp(1));
r=2*(M-1)/(P*M)*2^(B/(M-1));
else
r=1;
end
Pc=P*(1-r);
Pu=P*r/num_users;
Wc=zeros(M,num_users);

zf_F = (QH' * inv(QH * QH'));
Wc(:,1)=exp(i*2*pi*rand(4,1));
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
rate=rate1+rate2;
end
               
