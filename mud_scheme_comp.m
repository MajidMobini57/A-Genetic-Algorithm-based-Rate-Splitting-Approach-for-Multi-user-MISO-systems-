function [sum_rate, sum_rate_loss] = mud_scheme_comp(B,regions, Nt,K, P, BFType, CQIType,  MaxIter)
%K=number of users?

if strcmpi(BFType, 'P_ZF');
    BFType = 0; %perfect csit
    QuantizType=0;   %Perfect CSI
elseif strcmpi(BFType, 'ZF')
    BFType = 1; %RBF/PU2RC
     QuantizType=1;%Random vector quantization
else
    BFType = 2;  %Single-user beamforming
QuantizType=1; %Random vector quantization
end


if strcmpi(CQIType, 'Norm')
    CQIType = 0;
else
    CQIType = 1;    %SINR feedback 
end

    
if nargin < 8
    MaxIter = 100;
end
   

     temp_rate = 0;
    for its = 1:MaxIter

        H = 1/sqrt(2) * (randn(K, Nt) + 1i * randn(K, Nt)); %channel
        if BFType == 0 
             rate_temp = compute_ZF_rate(H,H,P);
        elseif BFType == 1  
            QuantH = quantiz_channels(H, regions, P, QuantizType, CQIType);
             rate_temp = compute_ZF_rate(QuantH,H, P);
        elseif BFType == 2 
            QuantH = quantiz_channels(H, regions, P, QuantizType, CQIType);
                          rate_temp = compute_P1_rate(QuantH,H,P,B,Nt);
%                           rate_temp = compute_P1_MU_rate(QuantH,H,P,B,Nt);


       
        end
 temp_rate=temp_rate+rate_temp;
    end; 
        sum_rate = temp_rate/MaxIter;
        
sum_rate_loss=0;
    end




%__EOF__
