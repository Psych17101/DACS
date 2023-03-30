function [FF,IFF] = PuckCriterion(sigma1,sigma2,sigma3,X_T,X_C,Y_T,Y_C,G12_t,nu12,E1,fiber)
FF =0;
IFF = 0;

if fiber == 'c'
    m = 1.1;
    % Incline parameters
    ptTll = 0.3;
    pcTll = 0.25;
    ptTT = 0.2;
    pcTT = 0.25;
    Ex = 230e9; % GPa  
elseif fiber == 'g'
    m = 1.3;
    % Incline parameters
    ptTll = 0.35;
    pcTll = 0.30;
    ptTT = 0.25;
    pcTT = 0.30;
    Ex = 75e9; % GPa
end

R_TTA = Y_C/(2*(1+pcTT));

if sigma1 > 0
    FF = 1/X_T*(sigma1 - (nu12-nu12*m*E1/Ex)*(sigma2+sigma3));
    fprintf("fiber failure Tension\n")
else 
    FF = 1/(-X_C)*(sigma1 - (nu12-nu12*m*E1/Ex)*(sigma2+sigma3));
    fprintf("Fibre failure Compression\n")
end

if sigma2 > 0  && Y_T > sigma2 %Condition Mode A
    R_Tll = G12_t;
    R_Tt = Y_T;
    A_ = (1/R_Tt - ptTll/R_Tll)*sigma2;
    B = sigma3/R_Tll;
    C = ptTll/R_Tll*sigma2;
    IFF = sqrt(A_^2 + B^2) + C;
    fprintf("Interfiber failure Mode A\n")
elseif sigma2 < 0 && -R_TTA < sigma2  %  Condition Mode B
    R_Tll = G12_t;
    A_ = (pcTll./R_Tll).*sigma2;
    B = sigma3./R_Tll;
    C = pcTll./R_Tll*sigma2;
    IFF = sqrt(A_.^2 + B.^2) + C;
    fprintf("Interfiber failure Mode B\n")
elseif sigma2 < -R_TTA && -Y_C <= sigma2 % Condition Mode C
    R_Tll = Y_C;
    A_ = sigma3./(2*(1+pcTT)*R_Tll);
    B = sigma2./R_Tll;
    C = R_Tll./-sigma2;
    IFF = (A_.^2 + B.^2)* C;
    thetafp = acosd(sqrt(1./2*(1+pcTT)*((R_TTA./R_Tll)*(sigma3./sigma2)+1)));
    fprintf("Interfiber failure Mode C\n")
    disp(thetafp)
else
    fprintf("invalid value of sigma2\n")
end



