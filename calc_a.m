function [a1 a2 a3 a4 a5] = calc_a(P, delta_a, delta_lambda, da_P, da_G, l, psy)
%CALC_A Summary of this function goes here
%   Detailed explanation goes here
    a1 = da_P'*da_P + psy*P'*P;
    a2 = 2*da_P'*(delta_a+da_G) + 2*psy*delta_lambda*P'*P;
    a3 = (delta_a+da_G)'*(delta_a+da_G)+psy*delta_lambda^2*P'*P-l^2;
    a4 = delta_a'*(delta_a+da_G);
    a5 = delta_a'*da_P;    

end

