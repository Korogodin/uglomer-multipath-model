%> @file test_CorrelFunc.m
%> @brief Корреляционный пик пробуем рисовать
%> @author Корогодин И.В.
%> @date   27 May 2011

globals;
nt = 1000;

dtm = -500:1:500;
TrueCorrFunc = ro(dtm);
RefCorrFunc = Amp_Ref(nt).*ro(dtm - Delta1(nt));
j = sqrt(-1);
absCorr = abs(TrueCorrFunc + RefCorrFunc*exp(j*2*pi*Delta1(nt)/lambda));

figure(908)
hA = gca;
plot(hA, dtm, TrueCorrFunc, dtm, RefCorrFunc, 'r');
hold(hA, 'on');
plot(hA, NaN, NaN, dtm, absCorr, 'LineWidth', 2)
hold(hA, 'off');



