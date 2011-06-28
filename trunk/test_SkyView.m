%> @file test_SkyView.m
%> @brief Пробуем функцию отрисовки SkyView
%> @author Корогодин И.В.
%> @date   24 May 2011

% Использует данные, оставшиеся от main-скрипта

globals;
nt = 3000;

lt = length(t);
Sky_x = nan(1, lt); Sky_y = nan(1, lt);
RefBeam_x = nan(1, lt); RefBeam_y = nan(1, lt);
for i = 1:lt
    % Координаты антенны в СК связанной с ЦЗ
    user.X = xa;
    user.Y = ya;
    user.Z = za + Re;

    if i == 10001
        disp('sdfsdf');
    end
    % Координаты спутника в СК связанной с ЦЗ
    st.X = xsv(i);
    st.Y = ysv(i);
    st.Z = zsv(i) + Re;

    ElAz = GetEaAz(user, st);
    Sky_x(i) = ElAz.X * 180 /pi;
    Sky_y(i) = ElAz.Y;
    if (zsv(i) <= za)
        Sky_y(i) = NaN;
        Sky_x(i) = NaN;
    end

    % Координаты точки отражения
    st.X = xo(i);
    st.Y = yo(i);
    st.Z = zo(i) + Re;
    ElAz = GetEaAz(user, st);
    RefBeam_x(i) = ElAz.X * 180 /pi;
    RefBeam_y(i) = ElAz.Y;
    if (zo(i) <= za)
        RefBeam_y(i) = NaN;
        RefBeam_x(i) = NaN;
    end
    
    
end
    Sky_y = Sky_y - pi*(Sky_x<0);  
    Alpha_sv = abs(Sky_x);
    Sky_x = 90 - abs(Sky_x); % ноль на внешнем радиусе, 90 в нуле   

    RefBeam_y = RefBeam_y - pi*(RefBeam_x<0);  
    Alpha_o = abs(RefBeam_x);
    RefBeam_x = 90 - abs(RefBeam_x); % ноль на внешнем радиусе, 90 в нуле   

% Экран
    st.X = xa + 5; st.Y = ya; st.Z = za + 2 + Re;
    ElAz = GetEaAz(user, st);
    Ekr_x(1) = ElAz.X*180/pi;  Ekr_y(1) = ElAz.Y;
    st.X = xa - 5;    st.Y = ya;    st.Z = za + 2 + Re;
    ElAz = GetEaAz(user, st);
    Ekr_x(2) = ElAz.X*180/pi;    Ekr_y(2) = ElAz.Y;
    Ekr_y = Ekr_y - pi*(Ekr_x<0);
    Ekr_x = 90 - abs(Ekr_x);


hhF = figure(904);
polar_my(gca, 0, 90, 'b');
hold on
polar_my(gca, Sky_y, Sky_x);
polar_my(gca, RefBeam_y, RefBeam_x, 'r');
polar_my(gca, Sky_y(nt), Sky_x(nt), '*');
polar_my(gca, RefBeam_y(nt), RefBeam_x(nt), 'r*');
polar_my(gca, Ekr_y, Ekr_x, 'k');
hold off