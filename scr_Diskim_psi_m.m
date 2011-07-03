%phi_m0 = -pi/2; 
%psi = pi/4; 

phi0 = rand(1,1)*2*pi;
psi_m = -pi/4;       

I1 = A_IQ*(cos(phi0 - psi) + k*cos(phi0 - phi_m0 - psi_m)) + n_mnoj*randn(1,1)*stdn_IQ;
Q1 = -A_IQ*(sin(phi0 - psi) + k*sin(phi0 - phi_m0 - psi_m)) + n_mnoj*randn(1,1)*stdn_IQ;
I2 = A_IQ*(cos(phi0 + psi) + k*cos(phi0 - phi_m0 + psi_m)) + n_mnoj*randn(1,1)*stdn_IQ; 
Q2 = -A_IQ*(sin(phi0 + psi) + k*sin(phi0 - phi_m0 + psi_m)) + n_mnoj*randn(1,1)*stdn_IQ;

Psi_m_extr = (-2:0.1:2)*pi;
s_Psi_m_extr = length(Psi_m_extr);
for j_psi_m_extr = 1:s_Psi_m_extr
    Ud(j_psi_m_extr) = Ud_psi_m_chief( I1, Q1, I2, Q2, -psi, Psi_m_extr(j_psi_m_extr), phi_m0, k ) / A_IQ / 2 * 40 ;
end

figure(997)
hold on
plot( (psi_m - Psi_m_extr), Ud)
ylabel('Ud_{\psi,m}');
hold off
grid on
drawnow
return