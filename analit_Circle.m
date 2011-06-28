%> @file analit_Circle.m
%> @brief Скрипт для визуализации аналитической задачи =)
%> @author Корогодин И.В.
%> @date   24 May 2011

% Использует данные, оставшиеся от main-скрипта

nt = 1; % Момент времени, для которого строим 

% Множество точек окружности:
pc = 0:0.01:2*pi;
x_circle = cos(pc)*sqrt(Rc2(nt)) + xa;
z_circle = sin(pc)*sqrt(Rc2(nt)) + za;

% Соответсвующие точки второго уравнения:
xo = x_circle;
zo = ( xo*(za-zsv(nt)) + xa*zsv(nt) - xsv(nt)*za ) / (xa - xsv(nt)) ;

% Аналитическое значение xo
xo_analit = ( xsv(nt) + xa*Rsva(nt)/Rao(nt) ) / (1 + Rsva(nt)/Rao(nt));

figure(991);
plot(x_circle, z_circle, xo, zo, [xo_analit xo_analit], [min(xo) max(xo)]);
xlabel('x, m')
ylabel('z, m')
grid on