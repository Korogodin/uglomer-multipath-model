function Ud = Ud_phi_m0_chief( I1, Q1, I2, Q2, psi, psi_m, phi_m0, k )
%UD_PSI_CHIEF Summary of this function goes here
%   Detailed explanation goes here

% I1 = Ib; I2 = Ia;
U1 = (I1 + I2)*cos(psi) + (Q2-Q1)*sin(psi) + ... % Xc
    k*( I1*cos(phi_m0+psi_m) + Q1*sin(phi_m0+psi_m) ) + ...
    k*( I2*cos(phi_m0-psi_m) + Q2*sin(phi_m0-psi_m) );

U2 = (Q1 + Q2)*cos(psi) + (I1-I2)*sin(psi) + ... % Xs
    k*( Q1*cos(phi_m0+psi_m) - I1*sin(phi_m0+psi_m) ) + ...
    k*( Q2*cos(phi_m0-psi_m) - I2*sin(phi_m0-psi_m) );

U1h = ...
    k*( -I1*sin(phi_m0+psi_m) + Q1*cos(phi_m0+psi_m) ) + ...
    k*( -I2*sin(phi_m0-psi_m) + Q2*cos(phi_m0-psi_m) );

U2h = ...
    k*( -Q1*sin(phi_m0+psi_m) - I1*cos(phi_m0+psi_m) ) + ...
    k*( -Q2*sin(phi_m0-psi_m) - I2*cos(phi_m0-psi_m) );

U = sqrt(U1^2 + U2^2);

Ud = U1*U1h + U2*U2h;

global A_IQ
Ud = Ud/U + k*2*A_IQ*sin(phi_m0)*cos(psi+psi_m);
end

