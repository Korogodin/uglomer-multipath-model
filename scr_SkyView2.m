Sky_y = Sky_y - pi*(Sky_x<0);  
Sky_x = 90 - abs(Sky_x); % ноль на внешнем радиусе, 90 в нуле   

RefBeam_y = RefBeam_y - pi*(RefBeam_x<0);  
RefBeam_x = 90 - abs(RefBeam_x); % ноль на внешнем радиусе, 90 в нуле   

% Экран
st.X = Screen_Width_l; st.Y = 0; st.Z = Re + 1.1*za;
ElAz = GetEaAz(user, st);
Ekr_x(1) = ElAz.X*180/pi;    Ekr_y(1) = ElAz.Y;
st.X = Screen_Width_l; st.Y = 0; st.Z = Screen_Hight + Re;
ElAz = GetEaAz(user, st);
Ekr_x(2) = ElAz.X*180/pi;  Ekr_y(2) = ElAz.Y;
st.X = -Screen_Width_r; st.Y = 0; st.Z = Screen_Hight + Re;
ElAz = GetEaAz(user, st);
Ekr_x(3) = ElAz.X*180/pi;    Ekr_y(3) = ElAz.Y;
st.X = -Screen_Width_r; st.Y = 0; st.Z = Re + 1.1*za;
ElAz = GetEaAz(user, st);
Ekr_x(4) = ElAz.X*180/pi;    Ekr_y(4) = ElAz.Y;

Ekr_y = Ekr_y - pi*(Ekr_x<0);
Ekr_x = 90 - abs(Ekr_x);