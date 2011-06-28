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