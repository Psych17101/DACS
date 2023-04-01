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
sigma_range = linspace(-Y_C, Y_T, 101); %N\m^2



% Definition of on distributed biaxial load
N_range = sigma_range*h; %N/m

% For each stress index in range increment.
for i = 1:length(sigma_range)
    for j = 1:length(sigma_range)
        % Biaxial stress
        sigma1 = sigma_range(i);
        sigma2 = sigma_range(j);
        F = [N_range(i);N_range(j);0];
        LR(i,j) = N_range(i)/N_range(j);
   
        % Calcuation of global strain for first ply failure
        strain_glo = A\F;
        % Calculation of global stresses
        stress_glo = Qbar*strain_glo; % global sigmaxx etc...
        %Initialisation of Failure criterion envelope
        maxstress_envelope2(i,j,:) = 0;
        puck_envelope2(i,j,:) = 0; 

        % Calculations of Failure Criterion for each ply 
        for l = 1:Nplies
            % Calculations of local Strains and Stresses
            [eps_loc] = strain_gtol(strain_glo,thetadb(l));
            [sigma_loc] = stress_gtol(stress_glo,thetadb(l));% ply i angle in radians, from bottom
            % Failure index with Maximum Stress criterion
            
            [FI_1(i,j,l),FI_2(i,j,l),FI_3(i,j,l)]= MaxStress(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t);
           

            % Failure index with Pucks Criterion
            [FF(i,j,l),IFF(i,j,l)] = PuckCriterion(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t,nu12,E1,'c','n');
            % Calculation of Ply failure
%             if FF(i,j,l) > 1 || IFF(i,j,l) >1
%                 puck_envelope2(i,j,l) = 1;
%             end
        end
        
            % Considering all the plies - find highest failure index of loading ration (i,j) 
            % Maximum Stress
            max_FI_1(i,j) = max(FI_1(i,j,:)); % Maximum FI_1 over all plies for loading ratio (i,j)
            max_FI_2(i,j) = max(FI_2(i,j,:));
            max_FI_3(i,j) = max(FI_3(i,j,:));
            % Highest Failure index according to maxmimum Stress
            max_FI_MS(i,j) = max([max_FI_1(i,j),max_FI_2(i,j),max_FI_3(i,j)]);
            
            % Puck Criterion
            max_FF(i,j) = max(FF(i,j,:));
            max_IFF(i,j) = max(IFF(i,j,:));
            % Highest Failure index according to Puck Criterion
            max_FI_PUCK(i,j) = max([max_FF(i,j),max_IFF(i,j)]);
       
            % Maximum Stress
            FPF(i,j,:) = F/max_FI_MS(i,j);
            Nx_FPF_MS(i,j) = FPF(i,j,1);
            Ny_FPF_MS(i,j) = FPF(i,j,2);
            
            % Puck Criterion
            FPF(i,j,:) = F/max_FI_PUCK(i,j);
            Nx_FPF_PUCK(i,j) = FPF(i,j,1);
            Ny_FPF_PUCK(i,j) = FPF(i,j,2);
         
        end
end


figure;
hold on
grid on
plot(Nx_FPF_MS(:,:),Ny_FPF_MS(:,:),'r*')
plot(Nx_FPF_PUCK(:,:),Ny_FPF_PUCK(:,:),'b*')
xlabel('N_x [N]');
ylabel('N_y [N]');
title('Biaxial Stress Failure Envelopes - Maximum Stress');
legend('Maximum Stress','Puck');


%% LPF
% Calculations of Failure index
i = 1;
% Initialisation of 1 iteration of Modulus
E_temp = E2;
iter = 0;
ply_index = 0;
ply_failure = zeros(1,Nplies);
F = [1;1;0];
temp1 = 0;

while E_temp ~= 0 && iter < 50 % While the failure index of our laminate is less than 1
    z = 0:h_ply:h;
    if temp1 < 10^(-25)
        F = F; % looped to creat different biaxial forces
    else
        F = F/temp1;
    end

    FI_LPF_1 = zeros(1,Nplies);
    FI_LPF_2 = zeros(1,Nplies);
    FI_LPF_3 = zeros(1,Nplies);

    % ABD reset
    A = zeros(3,3);
    B = zeros(3,3);
    D = zeros(3,3);
    Qbar = zeros(3,3);
    Sbar = zeros(3,3);

    % Calculation of Stiffness Matrix
    [S, Q] = ReducedComplianceStiffness(E1,E_temp,nu12,G12);
    [moduli]= [E1 E_temp nu12 G12];

    zbar = zeros(1, Nplies);

    for j = 1:length(ply_index)
        if ply_index ~= 0
            ply_failure(ply_index(j)) = 1;
        end
    end

    for l = 1:Nplies
        zbar(l) = -(h + h_ply)/2 + l*h_ply;
        % For each ply we calculate the ABD Matrix
        [Qbar(:,:,l),Sbar(:,:,l)] = QbarandSbar(thetadb(l),moduli);

        if ply_failure(l) == 1
            Q = ReducedStiffness(E_temp, nu12,nu21,G12);
            Qbar(:,:,l) = TransformedReducedQ(Q,thetadb(l));
        end

        A = A + Qbar(:,:,l) * (z(l+1)-z(l)) ; %N/m, right dimensions?
        B = B + (1/2) * Qbar(:,:,l) * (z(l+1)^2-z(l)^2); %N
        D = D + (1/3) * Qbar(:,:,l) * (z(l+1)^3-z(l)^3); %Nm
        ABD = [A B; A D];
    end
    
    A_test = A;
    invA_test = inv(A);

    % Calcuation of global strain for first ply failure
    strain_glo = inv(A)*F;
    % Calculation of global stresses
    stress_glo = Qbar(:,:,1)*strain_glo; % global sigmaxx etc.

    for l = 1:Nplies
        % Calculations of local Strains and Stresses
        [eps_loc] = strain_gtol(strain_glo,thetadb(l));
        [sigma_loc] = stress_gtol(stress_glo,thetadb(l));% ply i angle in radians, from bottom

        [FI_LPF_1(i,l),FI_LPF_2(i,l),FI_LPF_3(i,l)]= MaxStress(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t);
    end
    
    % Considering all the plies - find highest failure index of loading ration (i,j) 
    % Maximum Stress
    max_FI_1_LPF(i,j) = max(FI_LPF_1(i,j,:)); % Maximum FI_1 over all plies for loading ratio (i,j)
    max_FI_2_LPF(i,j) = max(FI_LPF_2(i,j,:));
    max_FI_3_LPF(i,j) = max(FI_LPF_3(i,j,:));
    % Highest Failure index according to maxmimum Stress
    max_FI_MS_LPF(i,j) = max([max_FI_1_LPF(i,j),max_FI_2_LPF(i,j),max_FI_3_LPF(i,j)]);

    if ply_failure == ones(1,Nplies)
        fprintf('LFP1 =')
        LPF1 = F;
        disp(LPF1)
        break
    end

    if max_FI_MS_LPF(i,j) == max_FI_1_LPF(i,j) % if fibre failure
        ply_index = find(FI_LPF_1(i,:) == FI_LPF_1(~ply_failure));
    elseif max_FI_MS_LPF(i,j) == max_FI_2_LPF(i,j)
        ply_index = find(FI_LPF_2(i,:) == FI_LPF_2(~ply_failure));
        E_temp = 0.1*E_temp;
    else
        ply_index = find(FI_LPF_3(i,:) == FI_LPF_3(~ply_failure));
        E_temp = 0.1*E_temp;
    end

    iter = iter + 1;
end




%% Functions
function [FF,IFF] = PuckCriterion(sigma1,sigma2,sigma3,X_T,X_C,Y_T,Y_C,G12_t,nu12,E1,fiber,print)
FF =0;
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


function y = ReducedStiffness(E2,NU12,NU21,G12)
%ReducedStiffness This function returns the reduced stiffness
% matrix for fiber-reinforced materials.
% There are four arguments representing four
% material constants. The size of the reduced
% stiffness matrix is 3 x 3.
y = [0, NU12*E2/(1-NU12*NU21), 0 ; NU12*E2/(1-NU12*NU21), E2/(1-NU12*NU21), 0 ; 0, 0, G12];
end

function y = TransformedReducedQ(Q,theta)
%Qbar This function returns the transformed reduced
% stiffness matrix "Qbar" given the reduced
% stiffness matrix Q and the orientation
% angle "theta".
% There are two arguments representing Q and "theta"
% The size of the matrix is 3 x 3.
% The angle "theta" must be given in degrees.
m = cosd(theta);
n = sind(theta);
T = [m*m n*n 2*m*n ; n*n m*m -2*m*n ; -m*n m*n m*m-n*n];
Tinv = [m*m n*n -2*m*n ; n*n m*m 2*m*n ; m*n -m*n m*m-n*n];
y = Tinv*Q*T;
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

function [Sbar, Qbar] = Tranformed_SCbar(S,Q,T,invT)
% This function returns the off axis reduced stiffness
% matrix for fiber-reinforced materials.
% There are four arguments representing four
% material constants. The size of the reduced
% stiffness and compliance matrix is 3 x 3.
Sbar = invT*S*T;
Qbar = invT*Q*T;
end


function [strain_glo] = strain_ltog(strain_loc,angle)
% Sine of the angle of the lamina
s=sind(angle);
% Cosine of the angle of the lamina
c=cosd(angle);
% Transformation matrix
T=[c^2, s^2, 2*s*c; s^2, c^2, -2*s*c; -s*c, s*c, c^2-s^2;];
% Calculating the global strain vector using the Tinv matrix
strain_glo=inv(T)*strain_loc;
% Calculating espxy
strain_glo(3)=strain_glo(3)*2;
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
end

function [stress_glo] = stress_ltog(stress_loc,angle)
% Sine of the angle of the lamina
s=sind(angle);
% Cosine of the angle of the lamina
c=cosd(angle);
% Transformation matrix
T=[c^2, s^2, 2*s*c; s^2, c^2, -2*s*c; -s*c, s*c, c^2-s^2;];
% Calculating the global strain vector using the Tinv matrix
stress_glo=inv(T)*stress_loc;
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

function [ps1,ps2,tmax,thetaps,thetass] = principalstresses(stress_glo)
% Stress in the x direction
sx=stress_glo(1);
% Stress in the y direction
sy=stress_glo(2);
% Stress in the x-y plane
txy=stress_glo(3);
% Longitudinal and Transverse Stresses
ps_1=((sx+sy)/2)+sqrt((((sx-sy)/2)^2)+((txy)^2));
ps_2=((sx+sy)/2)-sqrt((((sx-sy)/2)^2)+((txy)^2));
thetaps=atand((2*txy)/(sx-sy))/2;
if ps_1 > ps_2
    ps1=ps_1;
    ps2=ps_2;
else
    ps1=ps_2;
    ps2=ps_1;
end
% Shear Stresses
tmax=sqrt((((sx-sy)/2)^2)+(txy^2));
thetass=atand(-(sx-sy)/(2*txy))/2;
end
