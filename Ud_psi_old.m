function Ud = Ud_psi_old( I1, Q1, I2, Q2, psi )
%UD_PSI_OLD Summary of this function goes here
%   Detailed explanation goes here

Ud =   2*cos(psi)*(I1*Q2 - I2*Q1) - 2*sin(psi)*(I1*I2+Q1*Q2) ;

end

