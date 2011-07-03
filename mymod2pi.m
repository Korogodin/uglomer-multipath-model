%> @file mymod2pi.m
%> @brief mod to pm pi
%> @author Korogodin, I.V.
%> @date   24 May 2011
%> @todo 

function [ y ] = mymod2pi( x )
%MYMOD2PI Переводит число в интервал +-pi

y = mod(x+pi, 2*pi) - pi;
end

