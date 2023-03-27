function [FFT,FFC,IFF] = PuckCriterion(sigma1,sigma2,sigma3,X_T,X_C,nu12,E1,Ex,Y_T,Y_C,S,fiber)
FFT = 0; 
FFC = 0;
IFF = 0;

if fiber == 'c'
    m = 1.1;
    % Incline parameters
    ptTll = 0.3;
    pcTll = 0.25;
    ptTT = 0.2;
    pcTT = 0.25;
elseif fiber == 'g'
    m = 1.3;
    % Incline parameters
    ptTll = 0.35;
    pcTll = 0.30;
    ptTT = 0.25;
    pcTT = 0.30;
end

R_TTA = Y_C/(2*(1+pcTT));

if sigma1 > 0
    FFT = 1/X_T*(sigma1 - (nu12-nu12*m*Ex/E1)*(sigma2+sigma3));
else
    FFC = 1/(-X_C)*(sigma1 - (nu12-nu12*m*Ex/E1)*(sigma2+sigma3));
end

if sigma2 > 0  && Y_T > sigma2 %Condition Mode A
    R_Tll = S;
    R_Tt = Y_T;
    A = (1/R_Tt - ptTll/R_Tll)*sigma2;
    B = sigma3/R_Tll;
    C = ptTll/R_Tll*sigma2;
    IFF = sqrt(A^2 + B^2) + C;
    fprintf("Interfiber failure Mode A")
elseif sigma2 < 0 && -R_TTA < sigma2  %  Condition Mode B
    R_Tll = S;
    A = (pcTll/R_Tll)*sigma2;
    B = sigma3/R_Tll;
    C = pcTll/R_Tll*sigma2;
    IFF = sqrt(A^2 + B^2) + C;
    fprintf("Interfiber failure Mode A")
elseif sigma2 < -R_TTA && -Y_C < sigma2 % Condition Mode C
    R_Tll = Y_C;
    A = sigma3/(2*(1+pcTT)*R_Tll);
    B = sigma2/R_Tll;
    C = R_Tll/sigma2;
    IFF = (A^2 + B^2)* C;
    thetafp = acos(sqrt(1/2*(1+pcTT)*((R_TTA/R_Tll)*(sigma3/sigma2)+1)))
    fprintf("Interfiber failure Mode A")
    disp(thetafp)
else
    print("invalid value of sigma2")
end



