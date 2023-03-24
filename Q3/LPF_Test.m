clear; format;
format short g;

% Define material properties
E1   = 140.e9 ; % Pa - direction modulus of the lamina
E1_f = 80.e9; % Pa - directional modulus of fibres
nu12 = .3 ;
E2   = 10.e9  ; % Pa
E2_f = 20.e9; % Pa - tranvers modulus of fibres
G12  = 5.e9  ; % Pa
nu21 = nu12 * E2 / E1 ;
X_T = 1500.e6;
X_C = 1480.e6;
Y_T = 50.e6;
Y_C = 220.e6;
g_12t = 70.e6;
sigma_max = 200.e8;

%% Initialisation
% Laminate definition (plies of equal thickness)


% Define stress ranges for plotting
thetadt = [0 90 +45 -45 -45 +45 90 0]; % ply angles in degrees, from top??? Not bottom?
Nplies = length(thetadt);
thetadb = fliplr(thetadt); % ply angles in degrees, from bottom
h_ply  = 0.125*10^(-3); 
h      = Nplies * h_ply ;

sigma_range = linspace(-sigma_max, sigma_max, 101);
N_range = sigma_range*h/100;
max_FI = zeros(length(sigma_range));

% Calculations of Failure index
for i = 1:length(N_range)
        % Initialisation of 1 iteration of Modulus
        E_temp = E2;
        iter = 0;
        max_FI_1(i) = 0;
        index = 0;
        while E_temp ~= 0 && iter < 10 % While the failure index of our laminate is less than 1           
            z = 0:h_ply:h;
            if max_FI_1(i) == 0
                F = [N_range(i);0;0]; % looped to creat different biaxial forces
            else
                F = [N_range(i)/max_FI_1(i);0;0]; % looped to creat different biaxial forces
            end

            if index ~= 0
                thetadt(index) = [];
                thetadt(Nplies-index) = [];
                Nplies = length(thetadt);
                thetadb = fliplr(thetadt); % ply angles in degrees, from bottom
            end



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
            for l = 1:Nplies
                zbar(l) = -(h + h_ply)/2 + l*h_ply;
                % For each ply we calculate the ABD Matrix
                [Qbar,Sbar] = QbarandSbar(thetadb(l),moduli);
                A = A + Qbar * (z(l+1)-z(l)) ; %N/m, right dimensions?
                B = B + (1/2)*Qbar * (z(l+1)^2-z(l)^2); %N
                D = D + (1/3)*Qbar * (z(l+1)^3-z(l)^3); %Nm
                ABD = [A B; A D];
            end
            %disp(A)
            A_test = A;
            invA_test = inv(A);
            E_x = inv(h*invA_test(1,1));

            % Calcuation of global strain for first ply failure
            strain_glo = invA_test*F;
            max_fe(i) = 0;
            % Calculation of global stresses
            stress_glo = Qbar*strain_glo; % global sigmaxx etc..            

            clear FI_1
            clear FI_2
            clear FI_3
            clear FI



            for l = 1:Nplies
                % Calculations of local Strains and Stresses
                [eps_loc] = strain_gtol(strain_glo,thetadb(l));
                [sigma_loc] = stress_gtol(stress_glo,thetadb(l));% ply i angle in radians, from bottom
                sigma_loc = Q*eps_loc;

                % Failure index with Maximum Stress Criterion
                % For each pair of N(i) and N(j) there will be a maximum
                % failure criterion - from all the calculations of each failure
                % index, find the maximum one.
               

                [FI_1(i,l),FI_2(i,l),FI_3(i,l)]= MaxStress(sigma_loc,X_T,X_C,Y_T,Y_C,g_12t);
                FI = [FI_1(i,l),FI_2(i,l),FI_3(i,l)];

                % Failure index with Pucks Criterion
                fe(i,l) = fiberfailure(sigma_loc(1),sigma_loc(2),sigma_loc(3),'c',X_T,X_C,E_x,E1_f,nu12);  % (num_samples,layer)
                % Calculation of Ply failure
                if fe(i,l) > max_fe(i)
                    max_fe(i)= fe(i,l); %maximum failure index with relative forces
                    numberply(i,l) = l;
                end
            end
            % Considering all the plies - find highest index for each pair of
            max_FI_1(i) = max(FI_1(i,:)); % Strength Ratio
            max_FI_2(i) = max(FI_2(i,:)); % Strength Ratio

            % Creation of Forces (e.g coordinate on NxN) envelope for first ply failure
            F_FPF_MS1(i) = N_range(i)/max_FI_1(i);
            F_FPF_FF1(i) = N_range(i)/max_fe(i);


            if max_FI_1(i) > max_FI_2(i)
                index = find(max(FI_1(i,:)));
            else
                index = find(max(FI_2(i,:)));
            end

            if max_FI_1(i) < 1 && max_FI_2(i) < 1
                E_temp = 0.1*E_temp;
            else
                E_final(i) = E_temp;
                E_temp = 0;
            end

            iter = iter + 1;
            %disp(iter)
        end
end

%% Functions
function fe = fiberfailure(sigma1,sigma2,sigma3,fiber,X_T,X_C,Ex,E1,nu12)

    if fiber == 'c'
        m = 1.1;
    elseif fiber == 'g'
        m = 1.3;
    end

    if sigma1(1) > 0
        fe = 1/X_T*(sigma1 - (nu12-nu12*m*Ex/E1)*(sigma2+sigma3));
    else
        fe = 1/(-X_C)*(sigma1 - (nu12-nu12*m*Ex/E1)*(sigma2+sigma3));
    end

end

function [FI_1,FI_2,FI_3]= MaxStress(sigma,X_T,X_C,Y_T,Y_C,S_f)

if sigma(1)>=0
    FI_1 = sigma(1)/X_T;
else
    FI_1 = -sigma(1)/X_C;
end

if sigma(2)>=0
    FI_2 = sigma(2)/Y_T;
else
    FI_2 = -sigma(2)/Y_C;
end
FI_3 = abs(sigma(3))/S_f;
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
T=[c^2, s^2, 2*s*c; s^2, c^2, -2*s*c; -s*c, s*c, c^2-s^2;];
% local stress
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