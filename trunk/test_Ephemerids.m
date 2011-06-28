%> @file test_Ephemerids.m
%> @brief Тест функции расчета координат спутников
%> @author Корогодин И.В.
%> @date   24 May 2011


t = 0:1:(24*60*60);

X = nan(1,length(t));
Y = nan(1,length(t));
Z = nan(1,length(t));

Re = 6371000; %m

for i = 1:length(t)
    resi = ephemerids(1,t(i));
    X(i) = resi(1);
    Y(i) = resi(2);
    Z(i) = resi(3);
    if ~mod(i, fix(length(t)/100))
        fprintf('Done: %.0f %% \n', i/length(t)*100);
    end
end

hF = 0;
hF = figure(hF + 1);
plot3(X,Y,Z)
grid on

hF = figure(hF + 1);
plot(t,X)
grid on
xlabel('t, s');
ylabel('X, m');

hF = figure(hF + 1);
plot(t,Y)
grid on
xlabel('t, s');
ylabel('Y, m');

hF = figure(hF + 1);
plot(t,Z)
grid on
xlabel('t, s');
ylabel('Z, m');

hF = figure(hF + 1);
plot(t, ( sqrt(X.^2 + Y.^2 +Z.^2) - Re ) *1e-3 )
grid on
xlabel('t, s');
ylabel('h_{orb}, km');