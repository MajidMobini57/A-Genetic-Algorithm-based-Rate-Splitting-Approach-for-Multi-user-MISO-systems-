function rate = compute_ZF_rate(QH,H, P)


num_users = size(H, 1);
zf_F = (QH' * inv(QH * QH'));
for k = 1:num_users
    zf_F(:,k) = sqrt(P/num_users) * zf_F(:,k)/norm(zf_F(:,k));  
end;
gain_matrix = abs(H * zf_F).^2;
temp = 1+ sum(gain_matrix,2);
temp = diag(temp) * ones(num_users);
int_power = temp - gain_matrix;
SINR_matrix = gain_matrix ./ int_power;
rate = sum(log2(1 + diag(SINR_matrix)));
