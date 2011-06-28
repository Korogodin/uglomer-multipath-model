%> @file test_Freqmeter.m
%> @brief Оцениваем период ошибки многолучевости
%> @author Корогодин И.В.
%> @date   24 May 2011

% Использует данные, оставшиеся от main-скрипта

globals;

signErrPhi = ((ErrPhi > 0) - 0.5)*2;

ErrPeriod = nan(1,lt);
Peri = NaN; t1 = NaN;
SignOld = signErrPhi(1);
for i = 2:lt
    if ~isnan(ErrPhi(i))
        if SignOld ~= signErrPhi(i)
            Peri = 2*(t(i) - t1);
            t1 = t(i);
            SignOld = signErrPhi(i);
        end
    else
        Peri = NaN;
        t1 = NaN;
    end
    ErrPeriod(i) = Peri;
end

figure(907)
plot(t, ErrPeriod)
xlabel('t, s');
ylabel('MP Period, s');