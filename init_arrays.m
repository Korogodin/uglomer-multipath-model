Rc2 = nan(1, lt); % Квадрат радиуса точки отражения
Rsva = nan(1, lt); % Расстояние SV - антенна
Rao = nan(1, lt); % Расстояние Антенна - точка отражения
xsv = nan(1, lt); ysv = nan(1, lt); zsv = nan(1, lt); % Координаты SV
xo = nan(1, lt); zo = nan(1, lt); yo = zeros(1, lt); % Координаты ТО
cos_gamma = nan(1, lt);  Delta1 = nan(1, lt); Delta2 = nan(1, lt); %РХЛ
ErrPhi = nan(1, lt); % Ошибка в фазе
MP_flag = ones(1, lt); % Флаг наличия многолучевости
Sky_x = nan(1, lt); Sky_y = nan(1, lt); % SkyView прямой сигнал
RefBeam_x = nan(1, lt); RefBeam_y = nan(1, lt); % SkyView отраж. сигнал
Alpha_sv = nan(1, lt); Alpha_o = nan(1, lt); % Углы возвышения пр и отр
xpr = nan(1, lt); zpr = nan(1, lt); % Точка прямого сигнала на экране
true_signal_blocked_by_a_screen = nan(1, lt); 
sat_above_the_skyline = nan(1, lt);
direct_signal_is = nan(1, lt); % Прямой сигнал приходит на антенну
direct_signal_received = nan(1, lt); % Приходит на антенну выше нулевого угла места
ref_signal_is = nan(1, lt); % Аналогично с отраженным
ref_signal_received = nan(1, lt); % 
Amp_Ref = nan(1,lt); % Ослабление отраженного сигнала относительно прямого на входе антенны
Font_Size = 8;