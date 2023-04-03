clear; format;
close;
format short g;

%% Problems
% calculations of Puck are off - wrong equation?
% is E1_f correct? - related to puck
% Problem with plotting
% LPF should not be the same in negative or positive
% FPF should not be same in negative or positive
% FPF is should not be that small
% Why are most of the points inside the envelope


% Define material properties
E1   = 165.e9 ; % Pa - direction modulus of the lamina
E1_f = 20.e9; % Pa - directional modulus of fibres %% MAY BE WRONG
nu12 = .35 ;
E2   = 8.44*10^(9)  ; % Pa
E2_f = 20.e9; % Pa - tranvers modulus of fibres
G12  = 7.e9  ; % Pa
nu21 = nu12 * E2 / E1 ;
X_T = 1920.e6;
X_C = 1200.e6;
Y_T = 107.e6;
Y_C = 250.e6;
g_12t = 70.e6;
sigma_max = 1*10^4; % 2 MNm

%% Initialisation of Composites Model
% Laminate definition (plies of equal thickness)

thetadt = [0 90 +45 -45 -45 +45 90 0 0 90 +45 -45 -45 +45 90 0];% ply angles in degrees, from top??? Not bottom?
Nplies = length(thetadt);
thetadb = fliplr(thetadt); % ply angles in degrees, from bottom
h_ply  = 0.125*10^(-3);           % SI units, meters
h      = Nplies * h_ply ;

z = 0:h_ply:h;

% Creation of layering position in ply starting from bottom
for i = 1:Nplies
  zbar(i) = - (h + h_ply)/2 + i*h_ply;
end

% ABD Matrix Initialisation
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);


% Creation of Reduced Compliance and Stiffness Matrix
[S, Q] = ReducedComplianceStiffness(E1,E2,nu12,G12);
[moduli]= [E1 E2 nu12 G12];

% Creation of ABD Matrix 
for l = 1:Nplies
    % For each ply we calculate the ABD Matrix
[Qbar,Sbar(:,:,l)] = QbarandSbar(thetadb(l),moduli);
A = A + Qbar * (z(i+1)-z(i)) ; %N/m, right dimensions?
B = B + (1/2)* Qbar* (z(i+1)^2-z(i)^2); %N
D = D + (1/3)* Qbar * (z(i+1)^3-z(i)^3); %Nm
ABD = [A B; A D];
end

A_test = A;
invA = inv(A);

% Initialisation of Equivalent Modulus of Laminate
E_x = inv(h*invA(1,1));

%% Plot and Calculations biaxial stress failure envelopes for Puck and Max Stress criteria
% Define stress ranges for plotting
sigma_range = linspace(-2*Y_C, 3*Y_T, 101); %N\m^2
tau_range = linspace(0, 3*Y_T, 101); %N\m^2
% Definition of on distributed biaxial load
N2_range = sigma_range*h; %N/m
tau2_range = tau_range*h; %N/m

% For each stress index in range increment.
for i = 1:length(N2_range)
    for j = 1:length(tau2_range)
        % Biaxial stress
        F = [0;N2_range(i);tau2_range(j)];
        % Calcuation of global strain for first ply failure
        strain_glo = A\F;
        % Calculation of global stresses
        stress_glo = Qbar*strain_glo; % global sigmaxx etc...
        %Initialisation of Failure criterion envelope
        puck_envelope(i,j) = 0;
        % Calculations of Failure Criterion for each ply
        for l = 1:Nplies
            % Calculations of local Strains and Stresses
            [eps_loc] = strain_gtol(strain_glo,thetadb(l));
            [sigma_loc] = stress_gtol(stress_glo,thetadb(l));% ply i angle in radians, from bottom

            % Failure index with Pucks Criterion
            [FF(i,j,l),IFF(i,j,l)] = PuckCriterion(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t,nu12,E1,'c','n');
            % Calculation of Ply failure
            if FF(i,j,l) > 1 || IFF(i,j,l) >1
               puck_envelope(i,j) = 1;
            end
        end

    end
end


contour(sigma_range, tau_range, puck_envelope, [0 1], 'g')


%% Functions
function [FF,IFF] = PuckCriterion(sigma1,sigma2,sigma3,X_T,X_C,Y_T,Y_C,G12_t,nu12,E1,fiber,print)
FF = 0;
IFF = 0;

sigma3 = abs(sigma3);

% Determination of which fiber used
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
tau_12c = G12_t*sqrt(1+2*pcTT);

% Fiber failure (direction of the fibers)
if sigma1 > 0
    FF = 1/X_T*(sigma1 - (nu12-nu12*m*E1/Ex)*(sigma2+sigma3));
    if print == 'y'
        fprintf("Evaluating Fiber Tension\n")
    end
else 
    FF = 1/(-X_C)*(sigma1 - (nu12-nu12*m*E1/Ex)*(sigma2+sigma3));
    if print == 'y'
        fprintf("Evaluating Compression\n")
    end
end

% Matrix Failure
if sigma2 > 0 
    if sigma2 > Y_T || sigma3 > G12_t
        IFF = max(abs(sigma2/Y_T),sigma3/G12_t);
    elseif sigma2 < Y_T %Condition Mode A
    R_Tll = G12_t;
    R_Tt = Y_T;
    A_ = (1/R_Tt - ptTll/R_Tll)*sigma2;
    B = sigma3/R_Tll;
    C = ptTll/R_Tll*sigma2;
    IFF = sqrt(A_^2 + B^2) + C;
    if print == 'y'
        fprintf("Interfiber failure Mode A\n")
    end
    end
elseif sigma2 < -Y_C || sigma3 > tau_12c
    IFF =max(abs(sigma2/Y_C),sigma3/tau_12c);
elseif sigma2 <= 0 && -R_TTA < sigma2  %  Condition Mode B
    R_Tll = G12_t;
    A_ = (pcTll./R_Tll).*sigma2;
    B = sigma3./R_Tll;
    C = pcTll./R_Tll*sigma2;
    IFF = sqrt(A_.^2 + B.^2) + C;
    if print == 'y'
        fprintf("Interfiber failure Mode B\n")
    end
elseif sigma2 < -R_TTA && -Y_C <= sigma2 % Condition Mode C
    R_Tll = Y_C;
    A_ = sigma3./(2*(1+pcTT)*R_Tll);
    B = sigma2./R_Tll;
    C = R_Tll./-sigma2;
    IFF = (A_.^2 + B.^2)* C;
    thetafp = acosd(sqrt(1./2*(1+pcTT)*((R_TTA./R_Tll)*(sigma3./sigma2)+1)));
    if print == 'y'
        fprintf("Interfiber failure Mode C\n")
        disp(thetafp)
    end
else
    if print == 'y'
        fprintf("invalid value of sigma2\n")
    end
end

end


function [FI_1,FI_2,FI_3]= MaxStress(sigma1,sigma2,sigma3,X_T,X_C,Y_T,Y_C,S_f)

if sigma1>=0
    FI_1 = sigma1/X_T;
else
    FI_1 = -sigma1/X_C;
end

if sigma2>=0
    FI_2 = sigma2/Y_T;
else
    FI_2 = -sigma2/Y_C;
end
FI_3 = abs(sigma3)/S_f;
end


function [Qbar,Sbar] = QbarandSbar(angle,moduli)
% Sine of the angle of the lamina
s=sind(angle);
% Cosine of the angle of the lamina
c=cosd(angle);
% transformed compliance matrix values 
Sb11=((1/moduli(1))*c^4)+((2*(-moduli(3)/moduli(1))+(1/moduli(4)))*s^2*c^2)+((1/moduli(2))*s^4);
Sb12=((-moduli(3)/moduli(1))*(s^4+c^4))+(((1/moduli(1))+(1/moduli(2))-(1/moduli(4)))*s^2*c^2);
Sb16=((2*(1/moduli(1))-2*(-moduli(3)/moduli(1))-(1/moduli(4)))*s*c^3)-((2*(1/moduli(2))-2*(-moduli(3)/moduli(1))-(1/moduli(4)))*s^3*c);
Sb22=((1/moduli(1))*s^4)+((2*(-moduli(3)/moduli(1))+(1/moduli(4)))*s^2*c^2)+((1/moduli(2))*c^4);
Sb26=((2*(1/moduli(1))-2*(-moduli(3)/moduli(1))-(1/moduli(4)))*s^3*c)-((2*(1/moduli(2))-2*(-moduli(3)/moduli(1))-(1/moduli(4)))*s*c^3);
Sb66=(2*(2*(1/moduli(1))+2*(1/moduli(2))-4*(-moduli(3)/moduli(1))-(1/moduli(4)))*s^2*c^2)+((1/moduli(4))*(s^4+c^4));
% transformed compliance matrix and transformed reduced stiffness matrix
Sbar=[Sb11 Sb12 Sb16;Sb12 Sb22 Sb26;Sb16 Sb26 Sb66];
Qbar=inv(Sbar);
end


function [S, Q] = ReducedComplianceStiffness(E_1,E_2,v_12,G_12)
% This function returns the reduced stiffness
% matrix for fiber-reinforced materials.
% There are four arguments representing four
% material constants. The size of the reduced
% stiffness and compliance matrix is 3 x 3.
v_21 = v_12*E_2/E_1;
S = [1/E_1 -v_12/E_1 0 ;
    -v_12/E_1 1/E_2 0;
    0 0 1/G_12];
Q = [E_1/(1-v_12*v_21), v_12*E_2/(1-v_12*v_21), 0 ;
    v_12*E_2/(1-v_12*v_21), E_2/(1-v_12*v_21), 0 ; 
    0,0,G_12];

end

function [strain_loc] = strain_gtol(strain_glo,angle)
% Sine of the angle of the lamina
s=sind(angle);
% Cosine of the angle of the lamina
c=cosd(angle);
% Reuter Matrix
R=[1 0 0; 0  1  0;0  0  2];
% Transformation matrix
T=[c^2, s^2, 2*s*c; s^2, c^2, -2*s*c; -s*c, s*c, c^2-s^2;];
% local strain
strain_loc=R*T*inv(R)*strain_glo;
%strain_loc=T*strain_glo;

end

function [stress_loc] = stress_gtol(stress_glo,angle)
% Global stresses
% Sine of the angle of the lamina
s=sind(angle);
% Cosine of the angle of the lamina
c=cosd(angle);
% Transformation matrix
R = [1 0 0; 0  1  0;0  0  2];
% Transformation matrix
T = [c^2, s^2, 2*s*c; s^2, c^2, -2*s*c; -s*c, s*c, c^2-s^2;];
% local strain
stress_loc=T*stress_glo;
end
