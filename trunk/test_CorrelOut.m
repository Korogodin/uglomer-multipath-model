%> @file test_SkyView.m
%> @brief Пробуем отрисовывать отклики на выходе коррелятора
%> @author Корогодин И.В.
%> @date   24 May 2011

% Использует данные, оставшиеся от main-скрипта

for nt = 6000:7000;

complex_masht = 2; % Масштаб картинки
Amp_Ref = 0.5; % Амплитуда отраженного относительно прямого

% Единичный вертикальный вектор:
Vector_1_x = [0 0 0.015 0 -0.015]*complex_masht;
Vector_1_y = [0 1 0.85 1  0.85]*complex_masht;
Vector_1_comp = Vector_1_x + j*Vector_1_y;

% Прямой сигнал
Vector_True_Signal_x = Vector_1_x;
Vector_True_Signal_y = Vector_1_y;

% Окружность, описываемая вектором отраженного сигнала
pc = 0:0.05:2*pi;
Circle_Ref_Signal_x = Amp_Ref*cos(pc)*complex_masht;
Circle_Ref_Signal_y = (Amp_Ref*sin(pc) + 1)*complex_masht;

% Отраженный сигнал
lambda = 2.98e8 / 1.57542e9;
dPhi_Ref_Sig = 2*pi*Delta1(nt)/lambda; j = sqrt(-1);
Vector_Ref_Signal_x = real(Amp_Ref*Vector_1_comp*exp(j*dPhi_Ref_Sig));
Vector_Ref_Signal_y = imag(Amp_Ref*Vector_1_comp*exp(j*dPhi_Ref_Sig)) + 1*complex_masht;

% Суммарный сигнал
VectorSumAmp = abs(1 + Amp_Ref*exp(j*dPhi_Ref_Sig));
VectorSumAngle = angle(1 + Amp_Ref*exp(j*dPhi_Ref_Sig));
Vector_Sum_Signal_x = real(VectorSumAmp*Vector_1_comp*exp(j*VectorSumAngle));
Vector_Sum_Signal_y = imag(VectorSumAmp*Vector_1_comp*exp(j*VectorSumAngle));

figure(905)
plot(Vector_True_Signal_x, Vector_True_Signal_y, ... % Прямой сигнал
     Circle_Ref_Signal_x, Circle_Ref_Signal_y, '--r', ... % Окружность
     Vector_Ref_Signal_x, Vector_Ref_Signal_y, 'r', ... % Отраженный сигнал
     Vector_Sum_Signal_x, Vector_Sum_Signal_y); % Суммарный сигнал
xlim([-1.2 1.2]*complex_masht)
ylim([-0.2 2.2]*complex_masht)

pause(1)
end