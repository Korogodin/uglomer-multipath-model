Rc2 = nan(1, lt); Rc2_1 = nan(1, lt); Rc2_2 = nan(1, lt); % Квадрат радиуса точки отражения
Rsva = nan(1, lt); Rsva1 = nan(1, lt); Rsva2 = nan(1, lt); % Расстояние SV - антенна
Rao = nan(1, lt); Rao1 = nan(1, lt); Rao2 = nan(1, lt); % Расстояние Антенна - точка отражения
xsv = nan(1, lt); ysv = nan(1, lt); zsv = nan(1, lt); % Координаты SV
xo = nan(1, lt); zo = nan(1, lt); yo = zeros(1, lt); % Координаты ТО
xo1 = nan(1, lt); zo1 = nan(1, lt); yo1 = zeros(1, lt); % Координаты ТО для антенны 1
xo2 = nan(1, lt); zo2 = nan(1, lt); yo2 = zeros(1, lt); % Координаты ТО для антенны 2
cos_gamma = nan(1, lt);  
Delta1 = nan(1, lt); Delta1_1 = nan(1, lt);  Delta1_2 = nan(1, lt); % Разность хода лучей
Delta2 = nan(1, lt); %РХЛ
ErrPhi = nan(1, lt); ErrPhi1 = nan(1, lt); ErrPhi2 = nan(1, lt); % Ошибка в фазе
MP_flag = ones(1, lt); % Флаг наличия многолучевости
Sky_x = nan(1, lt); Sky_y = nan(1, lt); % SkyView прямой сигнал
RefBeam_x = nan(1, lt); RefBeam_y = nan(1, lt); % SkyView отраж. сигнал
Alpha_sv = nan(1, lt); Alpha_o = nan(1, lt); % Углы возвышения пр и отр
xpr = nan(1, lt); zpr = nan(1, lt); xpr1 = nan(1, lt); zpr1 = nan(1, lt); xpr2 = nan(1, lt); zpr2 = nan(1, lt); % Точка прямого сигнала на экране
true_signal_blocked_by_a_screen = nan(1, lt);  true_signal_blocked_by_a_screen1 = nan(1, lt);  true_signal_blocked_by_a_screen2 = nan(1, lt); 
sat_above_the_skyline = nan(1, lt); sat_above_the_skyline1 = nan(1, lt); sat_above_the_skyline2 = nan(1, lt);
direct_signal_is = nan(1, lt); direct_signal_is1 = nan(1, lt); direct_signal_is2 = nan(1, lt); % Прямой сигнал приходит на антенну
direct_signal_received = nan(1, lt); direct_signal_received1 = nan(1, lt); direct_signal_received2 = nan(1, lt); % Приходит на антенну выше нулевого угла места
ref_signal_is = nan(1, lt); ref_signal_is1 = nan(1, lt); ref_signal_is2 = nan(1, lt); % Аналогично с отраженным
ref_signal_received = nan(1, lt); ref_signal_received1 = nan(1, lt); ref_signal_received2 = nan(1, lt); % 
Amp_Ref = nan(1,lt); Amp_Ref1 = nan(1,lt); Amp_Ref2 = nan(1,lt); % Ослабление отраженного сигнала относительно прямого на входе антенны
Font_Size = 8;