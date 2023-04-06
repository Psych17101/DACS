% Question 2a

clear; format;
format short g;


%% Initialisation

Ex_2a = zeros(90,1);
Ey_2a = zeros(90,1);
nu_xy_2a = zeros(90,1);
nu_yx_2a = zeros(90,1);
Gxy_2a = zeros(90,1);

theta_in_2a = zeros(90,1);

Ex_flex_2a = zeros(90,1);
Ey_flex_2a = zeros(90,1);
nu_xy_flex_2a = zeros(90,1);
nu_yx_flex_2a = zeros(90,1);
Gxy_flex_2a = zeros(90,1);

for theta_in = 1:90
thetadt = [30 theta_in -theta_in 60 60 60 60 -theta_in theta_in 30]; % ply angles in degrees, from top??? Not bottom?
Nplies = length(thetadt);
thetadb = fliplr(thetadt); % ply angles in degrees, from bottom
h_ply  = 0.125*10^(-3);           % SI units, meters
h      = Nplies * h_ply ;

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

  A = A + Qbar * (z(i+1)-z(i)) ; %N/m,
  B = B + (1/2)*Qbar * (z(i+1)^2-z(i)^2); %N 
  D = D + (1/3)*Qbar * (z(i+1)^3-z(i)^3); %Nm

  NT = NT + Qbar * (z(i+1)-z(i));
  MT = MT + (1/2)*Qbar * (z(i+1)^2-z(i)^2) ;
  

  ABD = [A B; A D];
end

A_inv = inv(A);
B_inv = inv(B);
D_inv = inv(D);


% Question 2a calculations

theta_in_2a(theta_in,1) = theta_in;
Ex_2a(theta_in,1) = (A(1,1)*A(2,2) - A(1,2)^2)/(h*A(2,2));
Ey_2a(theta_in,1) = (A(1,1)*A(2,2) - A(1,2)^2)/(h*A(1,1));
nu_xy_2a(theta_in,1) = A(1,2)/A(2,2);
nu_yx_2a(theta_in,1) = A(1,2)/A(1,1);
Gxy_2a(theta_in,1) = A(3,3)/h;

Ex_flex_2a(theta_in,1) = 12/(h^3*D_inv(1,1));
Ey_flex_2a(theta_in,1) = 12/(h^3*D_inv(2,2));
nu_xy_flex_2a(theta_in,1) = -D_inv(1,2)/D_inv(1,1);
nu_yx_flex_2a(theta_in,1) = -D_inv(1,2)/D_inv(2,2);
Gxy_flex_2a(theta_in,1) = 12/(h^3*D_inv(3,3));



end

% Figures for Question 2a
figure(1)
plot(theta_in_2a, Ex_2a*10^-9, 'LineWidth',2)
xlabel('θ (degrees)')
ylabel('E_x (GPa)')
grid on
set(gca,'FontSize',14)
%[min_Ex, min_index] = min(Ex_2a);
%min_theta = theta_in_2a(min_index);
%hold on
%plot(min_theta, min_Ex*10^-9, 'ro', 'MarkerSize', 8, 'LineWidth', 2)
%annotation('textbox', [0.6, 0.15, 0.1, 0.1], 'String', {['min E_x = ' sprintf('%0.1f', min_Ex*10^-9) ' GPa'], ...
    %['at θ = ' num2str(min_theta) ' degrees']}, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FontSize', 12)

figure(2)
plot(theta_in_2a, Ey_2a*10^-9, 'LineWidth',2)
xlabel('θ (degrees)')
ylabel('E_y (GPa)')
grid on
set(gca,'FontSize',14)
%[min_Ey, min_index] = min(Ey_2a);
%min_theta = theta_in_2a(min_index);
%hold on
%plot(min_theta, min_Ey*10^-9, 'ro', 'MarkerSize', 8, 'LineWidth', 2)
%annotation('textbox', [0.4, 0.15, 0.1, 0.1], 'String', {['min E_y = ' sprintf('%0.1f', min_Ey*10^-9) ' GPa'], ...
    %['at θ = ' num2str(min_theta) ' degrees']}, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FontSize', 12)

figure(3)
plot(theta_in_2a, nu_xy_2a, 'LineWidth',2)
xlabel('θ (degrees)')
ylabel('ν_x_y')
grid on
set(gca,'FontSize',14)
%[max_nu_xy, max_index] = max(nu_xy_2a);
%max_theta = theta_in_2a(max_index);
%hold on
%plot(max_theta, max_nu_xy, 'ro', 'MarkerSize', 8, 'LineWidth', 2)
%annotation('textbox', [0.38, 0.75, 0.1, 0.1], 'String', {['max ν_x_y = ' sprintf('%0.3f', max_nu_xy)], ...
    %['at θ = ' num2str(max_theta) ' degrees']}, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FontSize', 12)



figure(4)
plot(theta_in_2a, nu_yx_2a, 'LineWidth',2)
xlabel('θ (degrees)')
ylabel('ν_y_x')
grid on
set(gca,'FontSize',14)
%[max_nu_yx, max_index] = max(nu_yx_2a);
%max_theta = theta_in_2a(max_index);
%hold on
%plot(max_theta, max_nu_yx, 'ro', 'MarkerSize', 8, 'LineWidth', 2)
%annotation('textbox', [0.6, 0.68, 0.1, 0.1], 'String', {['max ν_y_x = ' sprintf('%0.3f', max_nu_yx)], ...
    %['at θ = ' num2str(max_theta) ' degrees']}, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FontSize', 12)


figure(5)
plot(theta_in_2a, Gxy_2a*10^-9, 'LineWidth',2)
xlabel('θ (degrees)')
ylabel('G_x_y (GPa)')
grid on
set(gca,'FontSize',14)
%[max_Gxy, max_index] = max(Gxy_2a*10^-9);
%max_theta = theta_in_2a(max_index);
%hold on
%plot(max_theta, max_Gxy, 'ro', 'MarkerSize', 8, 'LineWidth', 2)
%annotation('textbox', [0.46, 0.75, 0.1, 0.1], 'String', {['max G_x_y = ' sprintf('%0.1f', max_Gxy) ' GPa'], ...
    %['at θ = ' num2str(max_theta) ' degrees']}, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FontSize', 12)


figure(6)
plot(theta_in_2a, Ex_flex_2a*10^-9, 'LineWidth',2)
xlabel('θ (degrees)')
ylabel('E_x (GPa)')
grid on
set(gca,'FontSize',14)

figure(7)
plot(theta_in_2a, Ey_flex_2a*10^-9, 'LineWidth',2)
xlabel('θ (degrees)')
ylabel('E_y (GPa)')
grid on
set(gca,'FontSize',14)

figure(8)
plot(theta_in_2a, nu_xy_flex_2a, 'LineWidth',2)
xlabel('θ (degrees)')
ylabel('ν_x_y')
grid on
set(gca,'FontSize',14)
%[max_nu_xy_flex, max_index] = max(nu_xy_flex_2a);
%max_theta = theta_in_2a(max_index);
%hold on
%plot(max_theta, max_nu_xy_flex, 'ro', 'MarkerSize', 8, 'LineWidth', 2)
%annotation('textbox', [0.36, 0.7, 0.1, 0.1], 'String', {['max ν_x_y = ' sprintf('%0.3f', max_nu_xy_flex)], ...
   % ['at θ = ' num2str(max_theta) ' degrees']}, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FontSize', 12)

figure(9)
plot(theta_in_2a, nu_yx_flex_2a, 'LineWidth',2)
xlabel('θ (degrees)')
ylabel('ν_y_x')
grid on
set(gca,'FontSize',14)
%[max_nu_yx_flex, max_index] = max(nu_yx_flex_2a);
%max_theta = theta_in_2a(max_index);
%hold on
%plot(max_theta, max_nu_yx_flex, 'ro', 'MarkerSize', 8, 'LineWidth', 2)
%annotation('textbox', [0.52, 0.7, 0.1, 0.1], 'String', {['max ν_y_x = ' sprintf('%0.3f', max_nu_yx_flex)], ...
    %['at θ = ' num2str(max_theta) ' degrees']}, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FontSize', 12)

figure(10)
plot(theta_in_2a, Gxy_flex_2a*10^-9, 'LineWidth',2)
xlabel('θ (degrees)')
ylabel('G_x_y (GPa)')
grid on
set(gca,'FontSize',14)
%[max_Gxy_flex, max_index] = max(Gxy_flex_2a*10^-9);
%max_theta = theta_in_2a(max_index);
%hold on
%plot(max_theta, max_Gxy_flex, 'ro', 'MarkerSize', 8, 'LineWidth', 2)
%annotation('textbox', [0.45, 0.7, 0.1, 0.1], 'String', {['max G_x_y = ' sprintf('%0.1f', max_Gxy_flex) ' GPa'], ...
   % ['at θ = ' num2str(max_theta) ' degrees']}, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FontSize', 12)





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