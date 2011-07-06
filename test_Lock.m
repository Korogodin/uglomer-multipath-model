globals;
clc
% set(handles.ed_Tlock, 'Enabled', 'off');
drawnow
% Tlock = fix(str2double(get(handles.ed_Tlock, 'String')));
Tlock = 1000;
if (Tlock < 2)
    Tlock = 30;
%     set(handles.ed_Tlock, 'String', num2str(Tlock));
end
if (nt + Tlock) > t(end)
    disp('Value of Tlock is too high');
%     set(handles.ed_Tlock, 'Enabled', 'on');
    return;
end

Tc = 0.005;
t_before_interp = nt:1:(nt+Tlock-1);
t_lock = nt:Tc:(nt+Tlock-1);
s_t_lock = length(t_lock);

Rsva1_int = interp1(t_before_interp, Rsva1(t_before_interp), t_lock);
Rsva2_int = interp1(t_before_interp, Rsva2(t_before_interp), t_lock);
Rsva0_int = (Rsva1_int + Rsva2_int)/2;
Phase0_int = (Rsva0_int(1) - Rsva0_int) / lambda * 2*pi;
OmegaPhase0_int = [diff(Phase0_int) diff(Phase0_int(end-1:end))] / Tc;
Psi_int = (Rsva1_int - Rsva2_int)/2  / lambda * 2*pi;

Delta1_1_int = interp1(t_before_interp, Delta1_1(t_before_interp), t_lock);
Delta1_2_int = interp1(t_before_interp, Delta1_2(t_before_interp), t_lock);
Delta1_0_int = (Delta1_1_int + Delta1_2_int)/2;

Phi_m0_int = Delta1_0_int / lambda * 2*pi;
Psi_m_int = - ( Delta1_2_int - Delta1_1_int )/2 / lambda * 2*pi;

stdn_IQ = 100;
qcno_dB = 45;
[A_IQ, qcno] = qcno_change(qcno_dB, stdn_IQ, Tc); 

F = [1 Tc Tc^2/2;
     0 1  Tc;
     0 0  1      ]; % Переходная матрица

Hf_psi = 1; % Hz, полоса
Kf_psi = nan(3,1); % Вектор-столбец коэффициентов фильтра
Kf_psi(3) = (1.2*Hf_psi)^3; % Коэффициенты непрерывной системы в установившемся режиме
Kf_psi(2) = 2*(Kf_psi(3))^(2/3);
Kf_psi(1) = 2*(Kf_psi(3))^(1/3);
Kf_psi = Kf_psi*Tc; % Переход к коэффициентам дискретной системы

Hf_phi_m0 = 0.2; % Hz, полоса
Kf_phi_m0 = nan(3,1); % Вектор-столбец коэффициентов фильтра
Kf_phi_m0(3) = (1.2*Hf_phi_m0)^3; % Коэффициенты непрерывной системы в установившемся режиме
Kf_phi_m0(2) = 2*(Kf_phi_m0(3))^(2/3);
Kf_phi_m0(1) = 2*(Kf_phi_m0(3))^(1/3);
Kf_phi_m0 = Kf_phi_m0*Tc; % Переход к коэффициентам дискретной системы

Hf_psi_m = 0.05; % Hz, полоса
Kf_psi_m = nan(3,1); % Вектор-столбец коэффициентов фильтра
Kf_psi_m(3) = (1.2*Hf_psi_m)^3; % Коэффициенты непрерывной системы в установившемся режиме
Kf_psi_m(2) = 2*(Kf_psi_m(3))^(2/3);
Kf_psi_m(1) = 2*(Kf_psi_m(3))^(1/3);
Kf_psi_m = Kf_psi_m*Tc; % Переход к коэффициентам дискретной системы

Xpsi_extr = [Psi_int(1); 0; 0];
Xphi_m0_extr = [Phi_m0_int(1); 0; 0];
Xpsi_m_extr = [Psi_m_int(1); 0; 0];

k = 0.2;
n_mnoj = 1;
psi_extr_j = nan(1, s_t_lock);  
phi_m0_extr_j = nan(1, s_t_lock); 
psi_m_extr_j = nan(1, s_t_lock); 
for j_t_lock = 1:s_t_lock
        
        phi0 = Phase0_int(j_t_lock);
        phi_m0 = Phi_m0_int(j_t_lock);
        psi_m = Psi_m_int(j_t_lock);
        psi = Psi_int(j_t_lock); 
        I1 = A_IQ*(cos(phi0 - psi) + k*cos(phi0 - phi_m0 - psi_m)) + n_mnoj*randn(1,1)*stdn_IQ;
        Q1 = -A_IQ*(sin(phi0 - psi) + k*sin(phi0 - phi_m0 - psi_m)) + n_mnoj*randn(1,1)*stdn_IQ;
        I2 = A_IQ*(cos(phi0 + psi) + k*cos(phi0 - phi_m0 + psi_m)) + n_mnoj*randn(1,1)*stdn_IQ; 
        Q2 = -A_IQ*(sin(phi0 + psi) + k*sin(phi0 - phi_m0 + psi_m)) + n_mnoj*randn(1,1)*stdn_IQ;
        
%         scr_Diskim_psi;
%         scr_Diskim_psi_m;
%         scr_Diskim_phi_m0;
%         return
        
        Ud_psi = -Ud_psi_chief( I1, Q1, I2, Q2, -Xpsi_extr(1), Xpsi_m_extr(1), Xphi_m0_extr(1), k ) / A_IQ / 2;
        Ud_psi_m = Ud_psi_m_chief( I1, Q1, I2, Q2, -Xpsi_extr(1),  Xpsi_m_extr(1), Xphi_m0_extr(1), k ) / A_IQ / 2 * 40 ;
        Ud_phi_m0 = Ud_phi_m0_chief( I1, Q1, I2, Q2, -Xpsi_extr(1), Xpsi_m_extr(1), Xphi_m0_extr(1), k ) / A_IQ * 12;
                
        Xpsi_extr = Xpsi_extr + Kf_psi*Ud_psi;
        psi_extr_j(j_t_lock) = Xpsi_extr(1);
        Xpsi_extr = F*Xpsi_extr;

        Xphi_m0_extr = Xphi_m0_extr + Kf_phi_m0*Ud_phi_m0;
        phi_m0_extr_j(j_t_lock) = Xphi_m0_extr(1);
        Xphi_m0_extr = F*Xphi_m0_extr;
        
        Xpsi_m_extr = Xpsi_m_extr + Kf_psi_m*Ud_psi_m;
        psi_m_extr_j(j_t_lock) = Xpsi_m_extr(1);
        Xpsi_m_extr = F*Xpsi_m_extr;        

        if ~mod(j_t_lock, fix(s_t_lock/10))
            fprintf('Lock: %.0f per cents \n', 100*j_t_lock/s_t_lock);
        end
end

hF = 20;
hF = figure(hF + 1);
plot(t_lock, rad2deg(psi_extr_j), t_lock, rad2deg(Psi_int));
ylabel('\Psi_{extr}, deg')
xlabel('Time, sec');

hF = figure(hF + 1);
ErrPsi = rad2deg(mymod2pi( 2*(psi_extr_j - Psi_int) ));
% hold on
plot(t_lock, ErrPsi);
hold off
ylabel('\delta \Psi, deg')
xlabel('Time, sec');
grid on
fprintf('Mean Error (2 sigma) of est Psi = %.2f deg \n', 2*sqrt(ErrPsi*ErrPsi'/length(ErrPsi)))

hF = figure(hF + 1);
plot(t_lock, rad2deg(phi_m0_extr_j), t_lock, rad2deg(Phi_m0_int));
ylabel('\phi_{m,0}, deg')
xlabel('Time, sec');
grid on

hF = figure(hF + 1);
ErrPhi_m0 = rad2deg(mymod2pi( 2*(phi_m0_extr_j - Phi_m0_int) ));
% hold on
plot(t_lock, ErrPhi_m0);
hold off
ylabel('\delta \phi_{m,0}, deg')
xlabel('Time, sec');
grid on
fprintf('Mean Error (2 sigma) of est Phi_m0 = %.2f deg \n', 2*sqrt(ErrPhi_m0*ErrPhi_m0'/length(ErrPhi_m0)))

hF = figure(hF + 1);
plot(t_lock, rad2deg(psi_m_extr_j), t_lock, rad2deg(Psi_m_int));
ylabel('\Psi_{m}, deg')
xlabel('Time, sec');
grid on

hF = figure(hF + 1);
ErrPsi_m = rad2deg(mymod2pi( 2*(psi_m_extr_j - Psi_m_int) ));
% hold on
plot(t_lock, ErrPsi_m);
hold off
ylabel('\delta \Psi_{m}, deg')
xlabel('Time, sec');
grid on
fprintf('Mean Error (2 sigma) of est Psi_m = %.2f deg \n', 2*sqrt(ErrPsi_m*ErrPsi_m'/length(ErrPsi_m)))
