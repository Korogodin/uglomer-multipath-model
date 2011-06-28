%> @file test_Earth.m
%> @brief Пробуем Земной шарик
%> @author Корогодин И.В.
%> @date   25 May 2011

% Использует данные, оставшиеся от main-скрипта

globals;

[xsp, ysp, zsp] = sphere;

figure(906)
% hp0 = mesh(Re*xsp, Re*ysp, Re*zsp);
hold on
hp0 = plot3(xsv/Re, ysv/Re, (zsv+Re)/Re);
hold off

hp1 = axesm ('globe','Grid', 'off');
view(60,60)
axis off

% Display a surface
load geoid
hp2 = meshm(geoid, geoidrefvec);

% Display coastline vectors
load coast
hp3 = plotm(lat, long);

% zdir = [0 0 1];
% for theta_earth = 0:360
%     drawnow
%     rotate(hp0, zdir, 1, [0 0 1])
%     rotate(hp1, zdir, 1, [0 0 1])
%     rotate(hp2, zdir, 1, [0 0 1])
%     rotate(hp3, zdir, 1, [0 0 1])
%     pause(0.1);
% end
