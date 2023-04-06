% Question 2b

clear; format;
format short g;


%% Initialisation

thetadt = [0 0 30 60 90 90 60 30 0 0]; % ply angles in degrees
Nplies = length(thetadt);
thetadb = fliplr(thetadt); % ply angles in degrees, from bottom
h_ply  = 0.125*10^(-3);           % SI units, meters
h      = Nplies * h_ply ;
epsilon_loc = zeros(3,Nplies);
sigma_loc_tot = zeros(3,Nplies);

z = -h/2:h_ply:h/2;

for i = 1:Nplies
  zbar(i) = - (h + h_ply)/2 + i*h_ply;
end

% Ply engineering properties examples
E1   = 140.e9 ; % Pa
nu12 = .3 ;
E2   = 10.e9  ; % Pa
G12  = 5.e9  ; % Pa
nu21 = nu12 * E2 / E1 ;
X_T = 1500.e6;
X_C = 1200.e6;
Y_T = 50.e6;
Y_C = 250.e6;
S_t = 70.e6;

[S, Q] = ReducedComplianceStiffness(E1,E2,nu12,G12);

% ABD
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

NT = zeros(3,1);
MT = zeros(3,1);
[moduli]= [E1 E2 nu12 G12];

for i = 1:Nplies
  % For each ply we calculate the ABD Matrix
  [Qbar,Sbar] = QbarandSbar(thetadb(i),moduli);

  A = A + Qbar * (z(i+1)-z(i)) ; %N/m
  B = B + (1/2)*Qbar * (z(i+1)^2-z(i)^2); %N 
  D = D + (1/3)*Qbar * (z(i+1)^3-z(i)^3); %Nm

  NT = NT + Qbar * (z(i+1)-z(i));
  MT = MT + (1/2)*Qbar * (z(i+1)^2-z(i)^2) ;
  
  ABD = [A B; A D];
end

A_inv = inv(A);
B_inv = inv(B);
D_inv = inv(D);
N_t= [220000;800;0];

eps_glo = A_inv*N_t;

z_middle = [-0.625 -0.5 -0.375 -0.25 -0.125 0.125 0.25 0.375 0.5 0.625];

    for i = 1:Nplies
        % For each ply we calculate the ABD Matrix
        theta  = thetadb(i)*pi/180;  % ply i angle in radians, from bottom
        [eps_loc] = strain_gtol(eps_glo,thetadb(i));
        sigma_loc = Q*eps_loc;
        epsilon_loc(:,i) = eps_loc;
        sigma_loc_tot(:,i) = sigma_loc;
    end

% Figures

figure(1)
plot(epsilon_loc(1,:),z_middle,'o-','LineWidth',2)
xlabel('ε_x')
ylabel('z (mm)')
ylim([-0.7, 0.7])
grid on;
set(gca,'FontSize',14)
hold on
line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineWidth', 2)
hold off



figure(2)
plot(epsilon_loc(2,:),z_middle,'o-','LineWidth',2)
xlabel('ε_y')
ylabel('z (mm)')
ylim([-0.7, 0.7])
grid on;
set(gca,'FontSize',14)
hold on
line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineWidth', 2)
hold off


figure(3)
plot(epsilon_loc(3,:),z_middle,'o-','LineWidth',2)
xlabel('ε_s')
ylabel('z (mm)')
ylim([-0.7, 0.7])
grid on;
set(gca,'FontSize',14)
hold on
line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineWidth', 2)
hold off

figure(4)
plot(sigma_loc_tot(1,:)*10^-6,z_middle,'o-','LineWidth',2)
xlabel('σ_x (MPa)')
ylabel('z (mm)')
ylim([-0.7, 0.7])
grid on;
set(gca,'FontSize',14)
hold on
line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineWidth', 2)
hold off

figure(5)
plot(sigma_loc_tot(2,:)*10^-6,z_middle,'o-','LineWidth',2)
xlabel('σ_y (MPa)')
ylabel('z (mm)')
ylim([-0.7, 0.7])
grid on;
set(gca,'FontSize',14)
hold on
line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineWidth', 2)
hold off

figure(6)
plot(sigma_loc_tot(3,:)*10^-6,z_middle,'o-','LineWidth',2)
xlabel('σ_s (MPa)')
ylabel('z (mm)')
ylim([-0.7, 0.7])
grid on;
set(gca,'FontSize',14)
hold on
line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineWidth', 2)
hold off


%% Functions

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

function [FI_1,FI_2,FI_3]= Failure1(sigma,X_T,X_C,Y_T,Y_C,S_f)

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