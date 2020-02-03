function QuantH = quantiz_channels(H, regions, P, QuantizType, CQIType)
% Quantizes channels H using a codebook of size "regions", using:
%   QuantizType quantization (0 - No Quantization, 1- Random Vector Quantization, 2 - Quantization Upper Bound, 3 - Scalar Quantization), and
%   scaled by the CQI value depending on CQIType (0 - Channel Norm, 1 - SINR computed using quantization error and SNR P)
% H and QuantH are of dimensions (Number of users) X (Number of base station antennas)

K = size(H, 1);
Nt = size(H, 2);
QuantH = zeros(K, Nt);
  %Choose quantized channels by generating random quantization errors
    temp = rand(K,1);

    if QuantizType == 0
        cos_sq = ones(1, K); 
    elseif QuantizType == 1 
        if log2(regions) <= 80%QUB
            cos_sq = 1 - (1-temp.^(1/regions)).^(1/(Nt-1));%baraye ma ine chon B=10....hamun h^ ast 
        else 
            cos_sq = double(vpa(1 - (1-temp.^(1/sym(regions))).^(1/(Nt-1))));
        end

    end
 
    for k = 1:K
        channel_gain = norm(H(k, :));
        rv = H(k, :)/channel_gain;
        temp = randn(Nt,1) + 1i * randn(Nt,1);
        temp = (temp - rv * temp * rv')';
        if CQIType == 0
            CQI = channel_gain^2;%Norm CQI
        else
            CQI = channel_gain^2*cos_sq(k) / (1 + P*channel_gain^2*(1-cos_sq(k))); %SINR CQI
        end
        QuantH(k,:) = sqrt(CQI)* (sqrt(cos_sq(k)) * rv + sqrt(1-cos_sq(k)) * temp/norm(temp));
    end

end %End if doscalar

