clc
globals;

Tlock = 100;
Tc = 0.005;
t_before_interp = nt:1:(nt+Tlock-1);
t_lock = nt:Tc:(nt+Tlock-1);
s_t_lock = length(t_lock);

Rsva1_int = interp1(t_before_interp, Rsva1(t_before_interp), t_lock);
Rsva2_int = interp1(t_before_interp, Rsva2(t_before_interp), t_lock);
Rsva0_int = (Rsva1_int + Rsva2_int)/2;
Phase0_int = (Rsva0_int(1) - Rsva0_int) / lambda * 2*pi;
OmegaPhase0_int = [diff(Phase0_int) diff(Phase0_int(end-1:end))] / Tc;
Psi21_int = (Rsva1_int - Rsva2_int)/2  / lambda * 2*pi;

Delta1_1_int = interp1(t_before_interp, Delta1_1(t_before_interp), t_lock);
Delta1_2_int = interp1(t_before_interp, Delta1_2(t_before_interp), t_lock);
Delta1_0_int = (Delta1_1_int + Delta1_2_int)/2;

Phi_m0_int = Delta1_0_int / lambda * 2*pi;
Psi_m_int = - ( Delta1_2_int - Delta1_1_int ) / lambda * 2*pi;

hF = 20;
hF = figure(hF+1);
plot(t_lock, Rsva0_int, t_lock, Rsva1_int, t_lock, Rsva2_int)
ylabel('Rsva')

hF = figure(hF+1);
plot(t_lock, Rsva1_int - Rsva2_int) % => psi_21!
ylabel('\Delta Rsva')

hF = figure(hF+1);
plot(t_lock, Phase0_int)
ylabel('Phase0')

hF = figure(hF+1);
plot(t_lock, OmegaPhase0_int)
ylabel('OmegaPhase0_int')

hF = figure(hF+1);
plot(t_lock, Psi21_int)
ylabel('Psi21_int')

hF = figure(hF+1);
plot(t_lock, Phi_m0_int)
ylabel('Phi_m0_int')

hF = figure(hF+1);
plot(t_lock, Psi_m_int)
ylabel('Psi_m_int')

hF = figure(hF+1);
plot(t_lock, Delta1_0_int, t_lock, Delta1_1_int, t_lock, Delta1_2_int);
ylabel('\Delta_{MP}')