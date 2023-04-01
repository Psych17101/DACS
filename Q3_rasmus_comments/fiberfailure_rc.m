function fe = fiberfailure(sigma1,sigma2,sigma3,fiber,X_T,X_C,E1_f,E1,v_12)

    if fiber == 'c'
        m = 1.1;
    elseif fiber == 'g'
        m = 1.3;
    end

    if sigma1 > 0
        fe = 1/X_T*(sigma1 - (v_12 - v_12*m*(E1/E1_f)*(sigma2 + sigma3)));
    elseif sigma1 <0
        fe = 1/X_C*(sigma1 - (v_12 - v_12*m*(E1/E1_f)*(sigma2 + sigma3)));
    end

end