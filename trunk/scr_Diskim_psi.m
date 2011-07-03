%phi_m0 = -pi/2; 
%psi_m = ;       

phi0 = rand(1,1)*2*pi;
psi = pi/4; 

I1 = A_IQ*(cos(phi0 - psi) + k*cos(phi0 - phi_m0 - psi_m)) + n_mnoj*randn(1,1)*stdn_IQ;
Q1 = -A_IQ*(sin(phi0 - psi) + k*sin(phi0 - phi_m0 - psi_m)) + n_mnoj*randn(1,1)*stdn_IQ;
I2 = A_IQ*(cos(phi0 + psi) + k*cos(phi0 - phi_m0 + psi_m)) + n_mnoj*randn(1,1)*stdn_IQ; 
Q2 = -A_IQ*(sin(phi0 + psi) + k*sin(phi0 - phi_m0 + psi_m)) + n_mnoj*randn(1,1)*stdn_IQ;

Psi_extr = (-2:0.1:2)*pi;
s_Psi_extr = length(Psi_extr);
for j_psi_extr = 1:s_Psi_extr
    Ud(j_psi_extr) = -Ud_psi_chief( I1, Q1, I2, Q2, -Psi_extr(j_psi_extr), psi_m, phi_m0, k ) / A_IQ / 2 ;
end

figure(996)
hold on
plot( (psi - Psi_extr), Ud)
hold off
grid on
ylabel('Ud_{\psi}');
drawnow
return