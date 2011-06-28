%> @file main.m
%> @brief Модель многолучевого распространения сигнала среднеорбитальной
%> спутниковой навигационной системы
%> @author Корогодин И.В.
%> @date   24 May 2011
%> @todo   Идет разработка...

minutes = 4*24*60;
t = 0:1:(minutes*60);
lt = length(t);

l = 20; h = 3;
xa = 0; ya = l; za = h;
Xa = [xa; ya; za];

Re = 6371000; %Средний радиус Земли по данными Википедии

Rc2 = nan(1, lt);
Rsva = nan(1, lt);
Rao = nan(1, lt);
xsv = nan(1, lt); ysv = nan(1, lt); zsv = nan(1, lt);
xo = nan(1, lt); zo = nan(1, lt); yo = zeros(1, lt);
cos_gamma = nan(1, lt);  Delta1 = nan(1, lt); Delta2 = nan(1, lt);
ErrPhi = nan(1, lt); % Ошибка в фазе
MP_flag = ones(1, lt); % Флаг наличия многолучевости
t0 = 0;
for i = 1:lt
    ephe_ce = ephemerids(9,t(i) + t0);
    ephe = ephe_ce - [0; 0; Re; 0; 0; 0]; % Переход из СК, связанной с центром Земли, в рабочую
    Xsv = ephe(1:3); xsv(i) = Xsv(1); ysv(i) = Xsv(2); zsv(i) = Xsv(3);
    Rsva(i) = norm(Xsv - Xa);
    Rc2(i) = ya^2*(Rsva(i)^2/ysv(i)^2 - 1);
    Rao(i) = abs(ya*Rsva(i)/ysv(i));
    xo(i) = ( xsv(i) + xa*Rsva(i)/Rao(i) ) / (1 + Rsva(i)/Rao(i));
    zo(i) = ( zsv(i) + za*Rsva(i)/Rao(i) ) / (1 + Rsva(i)/Rao(i));     Xo = [xo(i); 0; zo(i)];
    cos_gamma(i) = ((xo(i) - xa)*(xsv(i)-xa) + (0 - ya)*(ysv(i)-ya) + (zo(i) - za)*(zsv(i)-za) ) / Rao(i) / Rsva(i);
    Delta1(i) = Rao(i) * (1 - cos_gamma(i));
    Delta2(i) = norm(Xo - Xsv) + Rao(i) - Rsva(i);
    ErrPhi(i) = rad2deg( angle(1 + 0.5*exp(j*(2*pi*Delta1(i)/0.2))) );
    if (ysv(i) <= 0) || (zsv(i) <= 0)
        MP_flag(i) = 0;
        Delta1(i) = NaN;
        Delta2(i) = NaN;
        ErrPhi(i) = NaN;
        cos_gamma(i) = NaN;
        xo(i) = NaN;
        zo(i) = NaN;
        Rao(i) = NaN;
        Rc2(i) = NaN;
    end
    if ~mod(i, fix(length(t)/100))
        fprintf('Done: %.0f %% \n', i/lt*100);
    end
end

hF = 0;
hF = figure(hF + 1);
plot(t,Rsva)
grid on
xlabel('t, s');
ylabel('Rsva, m');

hF = figure(hF + 1);
plot(t,sqrt(Rc2))
grid on
xlabel('t, s');
ylabel('Rc, m');

hF = figure(hF + 1);
plot(t,Rao)
grid on
xlabel('t, s');
ylabel('Rao, m');

hF = figure(hF + 1);
plot(xo,zo, '.')
grid on
xlabel('xo, m');
ylabel('zo, m');

hF = figure(hF + 1);
plot(t,rad2deg(atan2(zo, xo)))
grid on
xlabel('t, s');
ylabel('angle_o, deg');

hF = figure(hF + 1);
plot(t, cos_gamma)
grid on
xlabel('t, s');
ylabel('cos(\gamma), deg');

hF = figure(hF + 1);
plot(t, Delta1, t, Delta2)
grid on
xlabel('t, s');
ylabel('\Delta, m');

hF = figure(hF + 1);
plot(t, ErrPhi)
grid on
xlabel('t, s');
ylabel('ErrPhi, deg');

hF = figure(hF + 1);
plot(t, xsv, t, ysv, t, zsv)
grid on
xlabel('t, s');
ylabel('x_{SV}, y_{SV}, z_{SV}, m');


