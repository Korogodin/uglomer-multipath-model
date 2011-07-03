% Блокируется ли луч экраном
if abs( Xsv(2) - Xa(2) ) > 1 
    xpr(i) = Xa(1) - Xa(2)*(Xsv(1) - Xa(1)) ./ ( Xsv(2) - Xa(2) );
    zpr(i) = Xa(3) - Xa(2)*(Xsv(3) - Xa(3)) ./ ( Xsv(2) - Xa(2) );
    if ( (xpr(i)<Screen_Width_l)&&(xpr(i)>(-Screen_Width_r))&& ...
            (zpr(i)>=0)&&(zpr(i)<=Screen_Hight) ) && (ysv(i) < 0)
        true_signal_blocked_by_a_screen(i) = 1;
    else
        true_signal_blocked_by_a_screen(i) = 0;
        xpr(i) = NaN; zpr(i) = NaN;
    end
else
    xpr(i) = NaN; zpr(i) = NaN;
    true_signal_blocked_by_a_screen(i) = 0;
end  


if abs( Xsv(2) - Xa1(2) ) > 1 
    xpr1(i) = Xa1(1) - Xa1(2)*(Xsv(1) - Xa1(1)) ./ ( Xsv(2) - Xa1(2) );
    zpr1(i) = Xa1(3) - Xa1(2)*(Xsv(3) - Xa1(3)) ./ ( Xsv(2) - Xa1(2) );
    if ( (xpr1(i)<Screen_Width_l)&&(xpr1(i)>(-Screen_Width_r))&& ...
            (zpr1(i)>=0)&&(zpr1(i)<=Screen_Hight) ) && (ysv(i) < 0)
        true_signal_blocked_by_a_screen1(i) = 1;
    else
        true_signal_blocked_by_a_screen1(i) = 0;
        xpr1(i) = NaN; zpr1(i) = NaN;
    end
else
    xpr1(i) = NaN; zpr1(i) = NaN;
    true_signal_blocked_by_a_screen1(i) = 0;
end  


if abs( Xsv(2) - Xa2(2) ) > 1 
    xpr2(i) = Xa2(1) - Xa2(2)*(Xsv(1) - Xa2(1)) ./ ( Xsv(2) - Xa2(2) );
    zpr2(i) = Xa2(3) - Xa2(2)*(Xsv(3) - Xa2(3)) ./ ( Xsv(2) - Xa2(2) );
    if ( (xpr2(i)<Screen_Width_l)&&(xpr2(i)>(-Screen_Width_r))&& ...
            (zpr2(i)>=0)&&(zpr2(i)<=Screen_Hight) ) && (ysv(i) < 0)
        true_signal_blocked_by_a_screen2(i) = 1;
    else
        true_signal_blocked_by_a_screen2(i) = 0;
        xpr2(i) = NaN; zpr2(i) = NaN;
    end
else
    xpr2(i) = NaN; zpr2(i) = NaN;
    true_signal_blocked_by_a_screen2(i) = 0;
end  