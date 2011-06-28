function varargout = fig_main(varargin)
% FIG_MAIN M-file for fig_main.fig
%      FIG_MAIN, by itself, creates a new FIG_MAIN or raises the existing
%      singleton*.
%
%      H = FIG_MAIN returns the handle to a new FIG_MAIN or the handle to
%      the existing singleton*.
%
%      FIG_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIG_MAIN.M with the given input arguments.
%
%      FIG_MAIN('Property','Value',...) creates a new FIG_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fig_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fig_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fig_main

% Last Modified by GUIDE v2.5 28-May-2011 18:39:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fig_main_OpeningFcn, ...
                   'gui_OutputFcn',  @fig_main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before fig_main is made visible.
function fig_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fig_main (see VARARGIN)

clc
global camz
camz = 0;
% Choose default command line output for fig_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fig_main wait for user response (see UIRESUME)
% uiwait(handles.fig_main);


% --- Outputs from this function are returned to the command line.
function varargout = fig_main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function ed_t0_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ed_t0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_h_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ed_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_l_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ed_l_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_Calc.
function pb_Calc_Callback(hObject, eventdata, handles)
scr_widg_off;
drawnow
globals;
read_edits(handles);

lt = length(t);

xa = 0; ya = l; za = h; % Координаты антенны в местной СК
Xa = [xa; ya; za]; 

Re = 6371000; %Средний радиус Земли по данными Википедии

    cos_nu = Re / (Re + h);
    h_xord = Re*cos_nu;
    r_xord = sqrt(Re^2 - h_xord^2);
    h_ost = Re + h - h_xord;
    tg_alpha = r_xord / h_ost; % минус tg угла линия горизонта - ось антенны
    Alpha_sv_0 = - (90 - rad2deg(atan(tg_alpha)));

init_arrays;
for i = 1:lt
    
    % Координаты SV
    ephe_ce = ephemerids(9,t(i) + t0);
    ephe = ephe_ce - [0; 0; Re; 0; 0; 0]; % Переход из СК, связанной с центром Земли, в рабочую
    Xsv = ephe(1:3); xsv(i) = Xsv(1); ysv(i) = Xsv(2); zsv(i) = Xsv(3);
    
    scr_true_signal_blocked_by_a_screen; % Блокируется ли луч экраном
    tg_alpha_sv = sqrt( (xsv(i) - xa).^2 + (ysv(i) - ya).^2 ) ./ (za - zsv(i));
    sat_above_the_skyline(i) = (tg_alpha_sv > tg_alpha) + (tg_alpha_sv <= 0); % Спутник выше горизонта
    direct_signal_is(i) = sat_above_the_skyline(i) * (~true_signal_blocked_by_a_screen(i));
    direct_signal_received(i) = direct_signal_is(i) * (zsv(i) > za);
    
    Rsva(i) = norm(Xsv - Xa);
    Rc2(i) = ya^2*(Rsva(i)^2/ysv(i)^2 - 1);
    Rao(i) = abs(ya*Rsva(i)/ysv(i));
    xo(i) = ( xsv(i) + xa*Rsva(i)/Rao(i) ) / (1 + Rsva(i)/Rao(i));
    zo(i) = ( zsv(i) + za*Rsva(i)/Rao(i) ) / (1 + Rsva(i)/Rao(i));     Xo = [xo(i); 0; zo(i)];
%     cos_gamma(i) = ((xo(i) - xa)*(xsv(i)-xa) + (0 - ya)*(ysv(i)-ya) + (zo(i) - za)*(zsv(i)-za) ) / Rao(i) / Rsva(i);
%     Delta2(i) = Rao(i) * (1 - cos_gamma(i));
    Delta1(i) = norm(Xo - Xsv) + Rao(i) - Rsva(i); j = sqrt(-1);
    Amp_Ref(i) = 0.75*ro(Delta1(i));
    ErrPhi(i) = rad2deg( angle(1 + Amp_Ref(i)*exp(j*(2*pi*Delta1(i)/0.2))) );
    
    if ( (xo(i)<Screen_Width_l)&&(xo(i)>(-Screen_Width_r))&& ...
            (zo(i)>=0)&&(zo(i)<=Screen_Hight) ) && (ysv(i) > 0) && ...
            direct_signal_is(i)
        ref_signal_is(i) = 1;
        if (zo(i) > za)
            ref_signal_received(i) = 1;
        else
            ref_signal_received(i) = 0;
        end
    else
        ref_signal_is(i) = 0;
        ref_signal_received(i) = 0;
        Delta1(i) = NaN;
        ErrPhi(i) = NaN;
        xo(i) = NaN;
        zo(i) = NaN;
        Rao(i) = NaN;
        Rc2(i) = NaN;
    end
    
    scr_SkyView;
    
    if ~mod(i, fix(length(t)/10))
        set(handles.pb_Calc, 'String', [ num2str(round(i/lt*100)) ' %']);
        set(handles.sl_time, 'Value', i);
        drawnow
    end
end

T1 = diff(Delta1);      size_T1 = length(T1);
T2 = diff(Delta1, 2);   size_T2 = length(T2);
T3 = diff(Delta1, 3);   size_T3 = length(T3);
T4 = diff(Delta1, 4);   size_T4 = length(T4);
T5 = diff(Delta1, 5);   size_T5 = length(T5);   
    
scr_SkyView2;

nt = lt;
plot_all(handles);

set(handles.sl_time, 'Value', nt);
set(handles.pb_Calc, 'String', 'Calc');
scr_widg_on;


function ed_Tmod_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ed_Tmod_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%> ======================================================================
%> @brief Функция чтения из edit-box'ов в локальные переменные
%> @param handles Структура указателей главной формы
% ======================================================================
function read_edits(handles)
globals;
Tmod = str2double(get(handles.ed_Tmod, 'String'));
if (Tmod <= 5)||(Tmod > 4*24*60*60)
    Tmod = 4*60*60;
    set(handles.ed_Tmod, 'String', num2str(Tmod));
end

t0 = str2double(get(handles.ed_t0, 'String')) - 1;

l = str2double(get(handles.ed_l, 'String'));
if (l<=0)||(l>5000)
    l = 20;
    set(handles.ed_l, 'String', num2str(l));
end

h = str2double(get(handles.ed_h, 'String'));
if (h<=0)||(h>1000)
    h = 3;
    set(handles.ed_h, 'String', num2str(h));
end

Screen_Width_l = str2double(get(handles.ed_Screen_width_left, 'String'));
if (Screen_Width_l<=0)||(Screen_Width_l>500)
    Screen_Width_l = 30;
    set(handles.ed_Screen_width_left, 'String', num2str(Screen_Width_l));
end

Screen_Width_r = str2double(get(handles.ed_Screen_width_right, 'String'));
if (Screen_Width_r<=0)||(Screen_Width_r>500)
    Screen_Width_r = 30;
    set(handles.ed_Screen_width_right, 'String', num2str(Screen_Width_r));
end

Screen_Hight = str2double(get(handles.ed_Screen_hight, 'String'));
if (Screen_Hight<=0)||(Screen_Hight>500)
    Screen_Hight = 30;
    set(handles.ed_Screen_hight, 'String', num2str(Screen_Hight));
end

t = 1:1:Tmod;
set(handles.sl_time, 'Max', max(t));
set(handles.sl_time, 'Value', min(t));
set(handles.sl_time, 'Min', min(t));



%> ======================================================================
%> @brief Функция отрисовки графика Rsva
%> @param handles Структура указателей главной формы, либо ноль, если
%> отрисовывать в отдельном окне
%> @param hF Флаг использования отдельного окна
% ======================================================================
function plot_axes_Rsva(handles, hF)
globals;
if hF 
    figure(hF_Rsva);
    hA = gca;
else
    hA = handles.axes_Rsva;
    set(hA, 'FontSize', Font_Size);
end
plot(hA, t, Rsva)
hold(hA, 'on'); plot(hA, t(nt), Rsva(nt), '*'); hold(hA, 'off');
grid(hA, 'on');
xlabel(hA, 't, s');
ylabel(hA, 'R_{SVA}, m');
title(hA, 'Distance SV-Antenna');


%> ======================================================================
%> @brief Функция отрисовки графика Rao
%> @param handles Структура указателей главной формы, либо ноль, если
%> отрисовывать в отдельном окне
%> @param hF Флаг использования отдельного окна
% ======================================================================
function plot_axes_Rao(handles, hF)
globals;
if hF 
    figure(hF_Rao);
    hA = gca;
else
    hA = handles.axes_Rao;
    set(hA, 'FontSize', Font_Size);        
end
plot(hA, t, Rao, 'r')
hold(hA, 'on'); plot(hA, t(nt), Rao(nt), 'r*'); hold(hA, 'off');
grid(hA, 'on');
xlabel(hA, 't, s');
ylabel(hA, 'R_{AO}, m');
title(hA, 'Distance Antenna-Reflection Point');


%> ======================================================================
%> @brief Функция отрисовки точки отражения на экране
%> @param handles Структура указателей главной формы, либо ноль, если
%> отрисовывать в отдельном окне
%> @param hF Флаг использования отдельного окна
% ======================================================================
function plot_axes_xozo(handles, hF)
globals;

% Очертание экрана
x_scr = [Screen_Width_l Screen_Width_l -Screen_Width_r -Screen_Width_r Screen_Width_l];
z_scr = [0 Screen_Hight Screen_Hight 0 0];

if hF 
    figure(hF_xozo);
    hA = gca;
else
    hA = handles.axes_xozo;
    set(hA, 'FontSize', Font_Size);    
end
plot(hA, xo, zo, 'r', x_scr, z_scr, 'k')
hold(hA, 'on'); 
plot(hA, xo(nt), zo(nt), 'r*');
plot(hA, xpr(nt), zpr(nt), 'b*'); 
hold(hA, 'off');
grid(hA, 'on');    
xlabel(hA, 'x_O, m');
ylabel(hA, 'z_O, m');
if hF
    legend(hA, 'Reflection point', 'Screen')
end
title(hA, 'Point on the screen');


%> ======================================================================
%> @brief Функция отрисовки графика разности хода лучей
%> @param handles Структура указателей главной формы, либо ноль, если
%> отрисовывать в отдельном окне
%> @param hF Флаг использования отдельного окна
% ======================================================================
function plot_axes_Delta1(handles, hF)
globals;
if hF 
    figure(hF_Delta1);
    hA = gca;
else
    hA = handles.axes_Delta1;
    set(hA, 'FontSize', Font_Size);
end
plot(hA, t, Delta1)
hold(hA, 'on'); plot(hA, nt, Delta1(nt), '*'); hold(hA, 'off');
grid(hA, 'on');    
xlabel(hA, 't, s');
ylabel(hA, '\Delta_{MP}, m');
title(hA, 'The difference of the beam''s paths');
if hF
    plot_Teylor(hA);
end

%> ======================================================================
%> @brief Функция отрисовки графика ошибки, вносимой в фазу многолучевостью
%> @param handles Структура указателей главной формы, либо ноль, если
%> отрисовывать в отдельном окне
%> @param hF Флаг использования отдельного окна
% ======================================================================
function plot_axes_ErrPhi(handles, hF)
globals;
if hF 
    figure(hF_ErrPhi);
    hA = gca;
else
    hA = handles.axes_ErrPhi;
    set(hA, 'FontSize', Font_Size);
end
plot(hA, t, ErrPhi); 
hold(hA, 'on');
plot(hA, nt, ErrPhi(nt), '*');
hold(hA, 'off');
grid(hA, 'on');    
xlabel(hA, 't, s');
ylabel(hA, '\Delta\phi_{MP}, deg');    
title(hA, 'Multipath phase error');


%> ======================================================================
%> @brief Функция отрисовки графика координат SV
%> @param handles Структура указателей главной формы, либо ноль, если
%> отрисовывать в отдельном окне
%> @param hF Флаг использования отдельного окна
% ======================================================================
function plot_axes_xsvysvzsv(handles, hF)
globals;
if hF 
    figure(hF_xsvysvzsv); hA = gca;
else
    hA = handles.axes_xsvysvzsv;
    set(hA, 'FontSize', Font_Size);
end
plot(hA, t, xsv, t, ysv, t, zsv)
hold(hA, 'on');
plot(hA, nt, xsv(nt), '*', nt, ysv(nt), '*', nt, zsv(nt), '*')
hold(hA, 'off');
if hF
    legend('x_{SV}', 'y_{SV}', 'z_{SV}');
end
title(hA, 'SpaceVehicle''s coordinates');
grid(hA, 'on');    
xlabel(hA, 't, s');
ylabel(hA, 'x_{SV}, y_{SV}, z_{SV}, m');  


%> ======================================================================
%> @brief Функция отрисовки комплексного сигнала
%> @param handles Структура указателей главной формы, либо ноль, если
%> отрисовывать в отдельном окне
%> @param hF Флаг использования отдельного окна
% ======================================================================
function plot_axes_Complex(handles, hF)
globals;
complex_masht = 1; % Масштаб картинки
j = sqrt(-1);

% Единичный вертикальный вектор:
Vector_1_x = [0 0 0.015 0 -0.015]*complex_masht;
Vector_1_y = [0 1 0.85 1  0.85]*complex_masht;
Vector_1_comp = Vector_1_x + j*Vector_1_y;

% Прямой сигнал
Vector_True_Signal_x = Vector_1_x;
Vector_True_Signal_y = Vector_1_y;

% Окружность, описываемая вектором отраженного сигнала
pc = 0:0.05:2*pi;
Circle_Ref_Signal_x = Amp_Ref(nt)*cos(pc)*complex_masht;
Circle_Ref_Signal_y = (Amp_Ref(nt)*sin(pc) + 1)*complex_masht;

% Отраженный сигнал
lambda = 2.98e8 / 1.57542e9;
dPhi_Ref_Sig = 2*pi*Delta1(nt)/lambda; 
Vector_Ref_Signal_x = real(Amp_Ref(nt)*Vector_1_comp*exp(j*dPhi_Ref_Sig));
Vector_Ref_Signal_y = imag(Amp_Ref(nt)*Vector_1_comp*exp(j*dPhi_Ref_Sig)) + 1*complex_masht;

% Суммарный сигнал
VectorSumAmp = abs(1 + Amp_Ref(nt)*exp(j*dPhi_Ref_Sig));
VectorSumAngle = angle(1 + Amp_Ref(nt)*exp(j*dPhi_Ref_Sig));
Vector_Sum_Signal_x = real(VectorSumAmp*Vector_1_comp*exp(j*VectorSumAngle));
Vector_Sum_Signal_y = imag(VectorSumAmp*Vector_1_comp*exp(j*VectorSumAngle));

if hF 
    figure(hF_Complex); hA = gca;
else
    hA = handles.axes_Complex;
    set(hA, 'FontSize', Font_Size);
end
plot(hA, Vector_True_Signal_x, Vector_True_Signal_y, ... % Прямой сигнал
     Circle_Ref_Signal_x, Circle_Ref_Signal_y, '--r', ... % Окружность
     Vector_Ref_Signal_x, Vector_Ref_Signal_y, 'r', ... % Отраженный сигнал
     Vector_Sum_Signal_x, Vector_Sum_Signal_y); % Суммарный сигнал
title(hA, 'Output signal of correlator');
xlim(hA, [-1.2 1.2]*complex_masht)
ylim(hA, [-0.2 2.2]*complex_masht)    
grid(hA, 'on');  


%> ======================================================================
%> @brief Функция отрисовки графика периода ошибки многолучевости
%> @param handles Структура указателей главной формы, либо ноль, если
%> отрисовывать в отдельном окне
%> @param hF Флаг использования отдельного окна
% ======================================================================
function plot_axes_Period(handles, hF)
globals;

signErrPhi = ((ErrPhi > 0) - 0.5)*2;

lt = length(t);
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

if hF 
    figure(hF_Period);
    hA = gca;
    hold(hA, 'on');
else
    hA = handles.axes_Period;
    set(hA, 'FontSize', Font_Size);
end
plot(hA, t, ErrPeriod); hold(hA, 'on'); 
plot(hA, t(nt), ErrPeriod(nt), '*'); 
hold(hA, 'off'); grid(hA, 'on');    
xlabel(hA,'t, s'); ylabel(hA,'MP Period, s');
title(hA, 'Multipath error period');


%> ======================================================================
%> @brief Функция отрисовки SkyView
%> @param handles Структура указателей главной формы, либо ноль, если
%> отрисовывать в отдельном окне
%> @param hF Флаг использования отдельного окна
% ======================================================================
function plot_axes_SkyView(handles, hF)
globals;

if hF 
    figure(hF_SkyView); hA = gca;
else
    hA = handles.axes_SkyView;
    set(hA, 'FontSize', Font_Size);
end
    polar_my(hA, 0, 90, 'b');
    hold(hA, 'on');
    polar_my(hA, Sky_y, Sky_x);
    polar_my(hA, RefBeam_y, RefBeam_x, 'r');
    polar_my(hA, Sky_y(nt), Sky_x(nt), '*');
    polar_my(hA, RefBeam_y(nt), RefBeam_x(nt), 'r*');
    polar_my(hA, Ekr_y, Ekr_x, 'k');
    hold(hA, 'off')
    grid(hA, 'on')
    title(hA, 'SkyView');



%> ======================================================================
%> @brief Функция отрисовки 3D "анимации": экран, антенна и т.д.
%> @param handles Структура указателей главной формы, либо ноль, если
%> отрисовывать в отдельном окне
%> @param hF Флаг использования отдельного окна
% ======================================================================
function plot_axes_3D(handles, hF)
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
if direct_signal_is(nt)
    true_signal_x = [xa xsv(nt)];
    true_signal_y = [ya ysv(nt)];
    true_signal_z = [za zsv(nt)];
else
    true_signal_x = NaN;
    true_signal_y = NaN;
    true_signal_z = NaN;
end    

% Падающий луч многолучевости
MP_inc_beam_x = [xo(nt) xsv(nt)];
MP_inc_beam_y = [0  ysv(nt)];
MP_inc_beam_z = [zo(nt) zsv(nt)];

% Отраженный луч многолучевости
MP_ref_beam_x = [xo(nt) xa];
MP_ref_beam_y = [0  ya];
MP_ref_beam_z = [zo(nt) za];

if ref_signal_is(nt)
    masht = 1.5*max([sqrt(Rc2(nt)), za, ya]);
else
    masht = 1.5*max([Screen_Width_l Screen_Width_r Screen_Hight, za, ya]);
end
    xpos = masht;    
    ypos = masht;    
    zpos = masht;
    xlim1 = -masht;
    xlim2 = masht;
    ylim1 = -0.1*masht;
    ylim2 = masht;
    zlim1 = -0.1*masht;
    zlim2 = masht;    

% Wall
wall_r = min([Screen_Width_r, abs(xlim1)]);
wall_l = min([Screen_Width_l, xlim2]);
wall_h = min([Screen_Hight, zlim2]);
wall_x = [-wall_r wall_l wall_l -wall_r -wall_r];
wall_y = [0 0 0 0 0];
wall_z = [0 0 wall_h wall_h 0];

if hF 
    hF = hF_3D;
    figure(hF); hA = gca;
else
    hA = handles.axes_3D;
    set(hA, 'FontSize', Font_Size);
end
    plot3(hA, x_circle, y_circle, z_circle, 'g', ...  % Окружность точки отражения
          graph_a_x, graph_a_y, graph_a_z, 'k', ... % Антенна
          true_signal_x, true_signal_y, true_signal_z, 'b', ... % Прямой сигнал
          wall_x, wall_y, wall_z, 'k', ... % Экран
          MP_inc_beam_x, MP_inc_beam_y, MP_inc_beam_z, 'r', ... % Падающий луч
          MP_ref_beam_x, MP_ref_beam_y, MP_ref_beam_z, 'r'); % Отраженный луч
    grid(hA, 'on');
    xlabel(hA, 'x');
    ylabel(hA, 'y');
    zlabel(hA, 'z');
    xlim(hA, [xlim1 xlim2]);
    ylim(hA, [ylim1 ylim2]);
    zlim(hA, [zlim1 zlim2]);
    camtarget(hA, [xa ya za]);
    campos(hA, [xpos, ypos, zpos]); 

    
%> ======================================================================
%> @brief Функция отрисовки графика угла места прямого и отраженного луча
%> @param handles Структура указателей главной формы, либо ноль, если
%> отрисовывать в отдельном окне
%> @param hF Флаг использования отдельного окна
% ======================================================================
function plot_axes_Angle(handles, hF)
globals;

if hF 
    figure(hF_Angle); hA = gca;
else
    hA = handles.axes_Angle;
    set(hA, 'FontSize', Font_Size);
end
plot(hA, t, Alpha_sv, t, Alpha_o, 'r', t, Alpha_sv_0*ones(1, length(t)), 'k')
hold(hA, 'on'); 
plot(hA, t(nt), Alpha_sv(nt), '*', t(nt), Alpha_o(nt), 'r*');
hold(hA, 'off');
grid(hA, 'on');   
xlabel(hA, 't, s');
ylabel(hA, '\alpha_e, deg');
if hF
    legend(hA, 'Direct', 'Reflected', 'Skyline')
end
title(hA, 'Angle of elevation');    

    
%> ======================================================================
%> @brief Функция отрисовки 3D графика Земли
%> @param handles Структура указателей главной формы, либо ноль, если
%> отрисовывать в отдельном окне
%> @param hF Флаг использования отдельного окна
% ======================================================================
function plot_axes_Earth(handles, hF)
globals;

if hF 
    figure(hF_Earth); 
    hA = gca;
    title(hA, 'Space View');
    hp0 = plot3(hA, xsv/Re, ysv/Re, (zsv+Re)/Re, 'c');
    hold(hA, 'on');
    hp01 = plot3(hA, xsv(nt)/Re, ysv(nt)/Re, (zsv(nt)+Re)/Re, 'b*');
    if direct_signal_is(nt)
        plot3(hA, [xsv(nt)/Re 0], [ysv(nt)/Re 0], [(zsv(nt)+Re)/Re 1+h/Re], '*-');
    end
    hold(hA, 'off');
    hp1 = axesm ('globe','Grid', 'off');
    %view(60,60)
    axis off
    %grid on
    % Display a surface
    load geoid
    hp2 = meshm(geoid, geoidrefvec);
    % Display coastline vectors
    load coast
    hp3 = plotm(lat, long);
    title(hA, 'Space View');
    camtarget(hA, [0 0 0]);
    campos(hA, [xsv(nt)/Re, ysv(nt)/Re, (zsv(nt)+Re)/Re]);    
    
else
    hA = handles.axes_Earth;
    cla(hA);
    set(handles.fig_main,'CurrentAxes', hA);
    set(hA, 'FontSize', Font_Size);
    hp0 = plot3(hA, xsv/Re, ysv/Re, (zsv+Re)/Re, 'c');
    hold(hA, 'on');
    hp01 = plot3(hA, xsv(nt)/Re, ysv(nt)/Re, (zsv(nt)+Re)/Re, 'b*');
    if direct_signal_is(nt)
        set(handles.fig_main,'CurrentAxes', hA);
        plot3(hA, [xsv(nt)/Re 0], [ysv(nt)/Re 0], [(zsv(nt)+Re)/Re 1+h/Re], '*-');
    end
    hold(hA, 'off');
    hp1 = axesm ('globe','Grid', 'off');
    %view(60,60)
    axis off
    %grid on
    % Display a surface
    load geoid
    hp2 = meshm(geoid, geoidrefvec);
    % Display coastline vectors
    load coast
    hp3 = plotm(lat, long);
    title(hA, 'Space View');
    camtarget(hA, [0 0 0]);
%     campos(hA, [xsv(nt)/Re, ysv(nt)/Re, (zsv(nt)+Re)/Re]);    
    campos(hA, [5, 5, 5]);    
    if camz == 0
        camzoom(hA, 3)
        camz = 1;
    end    
end
    
    
%> ======================================================================
%> @brief Функция отрисовки графика искаженной корреляционной функции
%> @param handles Структура указателей главной формы, либо ноль, если
%> отрисовывать в отдельном окне
%> @param hF Флаг использования отдельного окна
% ======================================================================
function plot_axes_Ro(handles, hF)
globals;

dtm = -500:1:500;
TrueCorrFunc = ro(dtm);
RefCorrFunc = Amp_Ref(nt).*ro(dtm - Delta1(nt));
j = sqrt(-1);
absCorr = abs(TrueCorrFunc + RefCorrFunc*exp(j*2*pi*Delta1(nt)/lambda));

if hF 
    figure(hF_Ro); 
    hA = gca;
else
    hA = handles.axes_Ro;
    set(hA, 'FontSize', Font_Size);
end
    plot(hA, dtm, TrueCorrFunc, dtm, RefCorrFunc, 'r');
    hold(hA, 'on'); 
    plot(hA, NaN, NaN, dtm, absCorr, 'LineWidth', 2)
    hold(hA, 'off');
    grid(hA, 'on'); 
    if hF 
        legend(hA, 'Direct', 'Reflected', 'Sum');
    end
    xlabel(hA, '\Delta\tau, s');
    ylabel(hA, '\rho(\Delta\tau)');
    title(hA, 'Correlation function');
    
    
% --- Executes on mouse press over axes background.
function axes_Rsva_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes on button press in pb_open_axes_Rsva.
function pb_open_axes_Rsva_Callback(hObject, eventdata, handles)
plot_axes_Rsva(0, 1);

% --- Executes on button press in pb_open_axes_Rao.
function pb_open_axes_Rao_Callback(hObject, eventdata, handles)
plot_axes_Rao(0, 1);

% --- Executes on button press in pb_open_axes_xozo.
function pb_open_axes_xozo_Callback(hObject, eventdata, handles)
plot_axes_xozo(0, 1);

% --- Executes on button press in pb_open_axes_Delta1.
function pb_open_axes_Delta1_Callback(hObject, eventdata, handles)
plot_axes_Delta1(0, 1);

% --- Executes on button press in pb_open_axes_ErrPhi.
function pb_open_axes_ErrPhi_Callback(hObject, eventdata, handles)
plot_axes_ErrPhi(0, 1);

% --- Executes on button press in pb_open_axes_xsvysvzsv.
function pb_open_axes_xsvysvzsv_Callback(hObject, eventdata, handles)
plot_axes_xsvysvzsv(0, 1);

% --- Executes on button press in pb_open_axes_Complex.
function pb_open_axes_Complex_Callback(hObject, eventdata, handles)
plot_axes_Complex(0, 1);

% --- Executes on button press in pb_open_axes_Period.
function pb_open_axes_Period_Callback(hObject, eventdata, handles)
plot_axes_Period(0, 1);

% --- Executes on button press in pb_open_axes_SkyView.
function pb_open_axes_SkyView_Callback(hObject, eventdata, handles)
plot_axes_SkyView(0, 1);

% --- Executes on button press in pb_open_axes_3D.
function pb_open_axes_3D_Callback(hObject, eventdata, handles)
plot_axes_3D(0, 1);

% --- Executes on button press in pb_open_axes_Angle.
function pb_open_axes_Angle_Callback(hObject, eventdata, handles)
plot_axes_Angle(0, 1);

% --- Executes on button press in pb_open_axes_Earth.
function pb_open_axes_Earth_Callback(hObject, eventdata, handles)
plot_axes_Earth(0, 1);

% --- Executes on button press in pb_open_axes_Ro.
function pb_open_axes_Ro_Callback(hObject, eventdata, handles)
plot_axes_Ro(0, 1);

% --- Executes on slider movement.
function sl_time_Callback(hObject, eventdata, handles)
globals;

nt = round(get(hObject, 'Value'));
plot_all(handles);

% --- Executes during object creation, after setting all properties.
function sl_time_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function plot_all(handles)
globals;
plot_axes_Rsva(handles, 0);
plot_axes_Rao(handles, 0);
plot_axes_xozo(handles, 0);
plot_axes_Delta1(handles, 0);
plot_axes_ErrPhi(handles, 0);
plot_axes_xsvysvzsv(handles, 0);
plot_axes_Complex(handles, 0);
plot_axes_Period(handles, 0);
plot_axes_SkyView(handles, 0);
plot_axes_3D(handles, 0);
plot_axes_Angle(handles, 0);
plot_axes_Earth(handles, 0);
plot_axes_Ro(handles, 0);
set(handles.txt_time, 'String', [num2str(nt) ' s']);


% --- Executes on button press in pb_PlayX1.
function pb_PlayX1_Callback(hObject, eventdata, handles)
globals;
    
if get(hObject, 'Value') == 1      
    scr_widg_off;
    set(handles.pb_PlayX1, 'Enable', 'on');
    drawnow
    for nt = nt:length(t)
        if get(hObject, 'Value') == 1
            pause(1);
            plot_all(handles);
            drawnow
        else
            break;
        end        
    end
    set(handles.sl_time, 'Value', nt);
    set(handles.pb_PlayX1, 'Value', 0);
end
scr_widg_on;


% --- Executes on button press in pb_PlayX10.
function pb_PlayX10_Callback(hObject, eventdata, handles)
globals;
    
if get(hObject, 'Value') == 1      
    scr_widg_off;
    set(handles.pb_PlayX10, 'Enable', 'on');
    drawnow
   
    for nt = nt:length(t)
        if get(hObject, 'Value') == 1
            pause(0.1);
            if ~mod(nt,10)
                plot_all(handles);
                drawnow
            end        
        else
            break;
        end        
    end
    set(handles.sl_time, 'Value', nt);
    set(handles.pb_PlayX10, 'Value', 0);
end
scr_widg_on;



% --- Executes on button press in pb_PlayX100.
function pb_PlayX100_Callback(hObject, eventdata, handles)
globals;
    
if get(hObject, 'Value') == 1      
    scr_widg_off;
    set(handles.pb_PlayX100, 'Enable', 'on');
    drawnow
    for nt = nt:100:length(t)
        if get(hObject, 'Value') == 1
            pause(1);
            plot_all(handles);
            drawnow
        else
            break;
        end        
    end
    set(handles.sl_time, 'Value', nt);
    set(handles.pb_PlayX100, 'Value', 0);
end
scr_widg_on;



function ed_Screen_hight_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ed_Screen_hight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_Screen_width_left_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ed_Screen_width_left_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_Screen_width_right_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ed_Screen_width_right_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_Teylor(hA)
globals;

T0nt = Delta1(nt);

if nt <= size_T1
    T1nt = T1(nt);
else
    T1nt = T1(nt-1);
end
if nt <= size_T2-1
    if nt >= 2
        T2nt = T2(nt-1);
    else
        T2nt = T2(nt);
    end
else
    T2nt = T2(nt-2);
end
if nt <= size_T3-2
    if nt >= 3
        T3nt = T3(nt-2);
    else
        T3nt = T3(nt);
    end
else
    T3nt = T3(nt-3);
end
if nt <= size_T4-3
    if nt >= 4
        T4nt = T4(nt-3);
    else
        T4nt = T4(nt);
    end
else
    T4nt = T4(nt-4);
end
if nt <= size_T5-4
    if nt >= 5
        T5nt = T5(nt-4);
    else
        T5nt = T5(nt);
    end
else
    T5nt = T5(nt-5);
end
c = 3e8;
if ~(isnan(T1nt)||isnan(T2nt)||isnan(T3nt)||isnan(T4nt)||isnan(T5nt))
text('position', [nt Delta1(nt)], 'fontsize', 8, 'string', ...
                      {[' Taylor series:'  ], ...
                       [' T1 = ' num2str(T1nt/c*1e12) ' ps/s'  ], ...
                       [' T2 = ' num2str(T2nt/c*1e12) ' ps/s^2'], ...
                       [' T3 = ' num2str(T3nt/c*1e12) ' ps/s^3'], ...
                       [' T4 = ' num2str(T4nt/c*1e12) ' ps/s^4'], ...
                       [' T5 = ' num2str(T5nt/c*1e12) ' ps/s^5']});
end


% --- Executes on button press in pb_PlayX1000.
function pb_PlayX1000_Callback(hObject, eventdata, handles)
globals;
    
if get(hObject, 'Value') == 1      
    scr_widg_off;
    set(handles.pb_PlayX1000, 'Enable', 'on');
    drawnow
    for nt = nt:1000:length(t)
        if get(hObject, 'Value') == 1
            pause(1);
            plot_all(handles);
            drawnow
        else
            break;
        end        
    end
    set(handles.sl_time, 'Value', nt);
    set(handles.pb_PlayX1000, 'Value', 0);
end
scr_widg_on;
