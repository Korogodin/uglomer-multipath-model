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
    

    if ( (xo1(i)<Screen_Width_l)&&(xo1(i)>(-Screen_Width_r))&& ...
            (zo1(i)>=0)&&(zo1(i)<=Screen_Hight) ) && (ysv(i) > 0) && ...
            direct_signal_is1(i)
        ref_signal_is1(i) = 1;
        if (zo1(i) > Xa1(3))
            ref_signal_received1(i) = 1;
        else
            ref_signal_received1(i) = 0;
        end
    else
        ref_signal_is1(i) = 0;
        ref_signal_received1(i) = 0;
        Delta1_1(i) = NaN;
        ErrPhi1(i) = NaN;
        xo1(i) = NaN;
        zo1(i) = NaN;
        Rao1(i) = NaN;
        Rc2_1(i) = NaN;
    end    
    

    if ( (xo2(i)<Screen_Width_l)&&(xo2(i)>(-Screen_Width_r))&& ...
            (zo2(i)>=0)&&(zo2(i)<=Screen_Hight) ) && (ysv(i) > 0) && ...
            direct_signal_is2(i)
        ref_signal_is2(i) = 1;
        if (zo2(i) > Xa2(3))
            ref_signal_received2(i) = 1;
        else
            ref_signal_received2(i) = 0;
        end
    else
        ref_signal_is2(i) = 0;
        ref_signal_received2(i) = 0;
        Delta1_2(i) = NaN;
        ErrPhi2(i) = NaN;
        xo2(i) = NaN;
        zo2(i) = NaN;
        Rao2(i) = NaN;
        Rc2_2(i) = NaN;
    end    