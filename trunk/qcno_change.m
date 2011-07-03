%> @file qcno_change.m
%> @brief qcno_change function: SNR changing
%> @author Korogodin, I.V.
%> @date   24 May 2011
%> @todo 

function [A_IQ qcno] =  qcno_change(qcno_dB, stdn_IQ, Tc)
%QCNO_CHANGE Функция сопоставляет массиву qcno_dB (dBHz) -> 
% массив qcno, массив A_IQ

qcno = 10.^(qcno_dB/10);
A_IQ = stdn_IQ .* sqrt(2 * qcno * Tc);

end

