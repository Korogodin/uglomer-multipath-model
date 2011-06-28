% Расчет SkyView

% Координаты антенны в СК связанной с ЦЗ
user.X = xa;
user.Y = ya;
user.Z = za + Re;

% Координаты спутника в СК связанной с ЦЗ
st.X = xsv(i);
st.Y = ysv(i);
st.Z = zsv(i) + Re;

ElAz = GetEaAz(user, st);
Sky_x(i) = ElAz.X * 180 /pi;
Sky_y(i) = ElAz.Y;
Alpha_sv(i) = Sky_x(i);    
if ~direct_signal_received(i)
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
Alpha_o(i) = RefBeam_x(i);
if ~ref_signal_received(i)
    RefBeam_y(i) = NaN;
    RefBeam_x(i) = NaN;
end    
