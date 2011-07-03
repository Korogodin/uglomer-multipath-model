%psi = pi/4; 
%psi_m = pi/4;   
phi0 = rand(1,1)*2*pi;

phi_m0 = pi/3; 

I1 = A_IQ*(cos(phi0 - psi) + k*cos(phi0 - phi_m0 - psi_m)) + n_mnoj*randn(1,1)*stdn_IQ;
Q1 = -A_IQ*(sin(phi0 - psi) + k*sin(phi0 - phi_m0 - psi_m)) + n_mnoj*randn(1,1)*stdn_IQ;
I2 = A_IQ*(cos(phi0 + psi) + k*cos(phi0 - phi_m0 + psi_m)) + n_mnoj*randn(1,1)*stdn_IQ; 
Q2 = -A_IQ*(sin(phi0 + psi) + k*sin(phi0 - phi_m0 + psi_m)) + n_mnoj*randn(1,1)*stdn_IQ;

Phi_m0_extr = (-2:0.1:2)*pi;
s_Phi_m0_extr = length(Phi_m0_extr);
for j_phi_m0_extr = 1:s_Phi_m0_extr
    Ud(j_phi_m0_extr) = Ud_phi_m0_chief( I1, Q1, I2, Q2, -psi, psi_m, Phi_m0_extr(j_phi_m0_extr), k ) / A_IQ * 12 ;
end

figure(998)
hold on
plot( (phi_m0 - Phi_m0_extr), Ud)
ylabel('Ud_{\phi,m0}');
hold off
grid on
drawnow
return