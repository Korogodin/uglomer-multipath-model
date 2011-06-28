%> @file test_3D.m
%> @brief Пробуем функцию отрисовки 3D картинки
%> @author Корогодин И.В.
%> @date   24 May 2011

% Использует данные, оставшиеся от main-скрипта

nt = 1000;
globals;

% Множество точек окружности:
pc = 0:0.05:2*pi;
x_circle = cos(pc)*sqrt(Rc2(nt)) + xa;
z_circle = sin(pc)*sqrt(Rc2(nt)) + za; z_circle = z_circle.*(z_circle>0);
y_circle = zeros(1, length(pc));


% Антенна
graph_a_x = [0 0 0.2 -0.2 0 0 0 0];
graph_a_y = [ya ya ya ya ya ya+0.2 ya-0.2 ya];
graph_a_z = [0 za za+0.2 za+0.2 za za+0.2 za+0.2 za];

% Прямой сигнал
true_signal_x = [xa xsv(nt)];
true_signal_y = [ya ysv(nt)];
true_signal_z = [za zsv(nt)];

% Wall
masht = 1.5*sqrt(Rc2(nt)) + za;
wall_w = 1.3*sqrt(Rc2(nt)) + za;
wall_h = 1.2*sqrt(Rc2(nt)) + za;
wall_x = [-wall_w wall_w wall_w -wall_w -wall_w];
wall_y = [0 0 0 0 0];
wall_z = [0 0 wall_h wall_h 0];

% Падающий луч многолучевости
MP_inc_beam_x = [xo(nt) xsv(nt)];
MP_inc_beam_y = [0  ysv(nt)];
MP_inc_beam_z = [zo(nt) zsv(nt)];

% Отраженный луч многолучевости
MP_ref_beam_x = [xo(nt) xa];
MP_ref_beam_y = [0  ya];
MP_ref_beam_z = [zo(nt) za];

% Нормаль
normal_x = [xo(nt) xo(nt)];
normal_y = [0 0+0.4];
normal_z = [zo(nt) zo(nt)];


figure(902)
plot3(x_circle, y_circle, z_circle, 'g', ...  % Окружность точки отражения
      graph_a_x, graph_a_y, graph_a_z, 'k', ... % Антенна
      true_signal_x, true_signal_y, true_signal_z, 'b', ... % Прямой сигнал
      wall_x, wall_y, wall_z, 'k', ... % Экран
      MP_inc_beam_x, MP_inc_beam_y, MP_inc_beam_z, 'r', ... % Падающий луч
      MP_ref_beam_x, MP_ref_beam_y, MP_ref_beam_z, 'r'); % Отраженный луч
grid on
xlabel('xlabel')
ylabel('ylabel')
zlabel('zlabel')
xlim([-masht masht]);
ylim([-0.1*masht masht]);
zlim([-0.1*masht masht]);
campos([3*wall_w,2*ya,1.5*wall_h]);
% zdir = [0 0 1];
% rotate(902,zdir,0,Xa')

figure(903)
plot3(x_circle, y_circle, z_circle, 'g', ...  % Окружность точки отражения
      MP_inc_beam_x, MP_inc_beam_y, MP_inc_beam_z, 'r', ... % Падающий луч
      MP_ref_beam_x, MP_ref_beam_y, MP_ref_beam_z, 'r', ... % Отраженный луч
      normal_x, normal_y, normal_z, 'k'); % Нормаль 
grid on
xlabel('xlabel')
ylabel('ylabel')
zlabel('zlabel')
xlim([xo(nt)-1 xo(nt)+1]);
ylim([0-1 0+1]);
zlim([zo(nt)-1 zo(nt)+1]);
%campos([3*wall_w,2*ya,1.5*wall_h]);