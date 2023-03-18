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
sigma_max = 50.e8;

%% Initialisation
% Laminate definition (plies of equal thickness)

thetadt = [0 90 +45 -45 -45 +45 90 0]; % ply angles in degrees, from top??? Not bottom?
Nplies = length(thetadt);
thetadb = fliplr(thetadt); % ply angles in degrees, from bottom
h_ply  = 0.125*10^(-3);           % SI units, meters
h      = Nplies * h_ply ;

z = 0:h_ply:h;

for i = 1:Nplies
  zbar(i) = - (h + h_ply)/2 + i*h_ply;
end
 
% ABD
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

[S, Q] = ReducedComplianceStiffness(E1,E2,nu12,G12);
[moduli]= [E1 E2 nu12 G12];

for l = 1:Nplies
    % For each ply we calculate the ABD Matrix
[Qbar,Sbar] = QbarandSbar(thetadb(l),moduli);
A = A + Qbar * (z(i+1)-z(i)) ; %N/m, right dimensions?
B = B + (1/2)*Qbar * (z(i+1)^2-z(i)^2); %N
D = D + (1/3)*Qbar * (z(i+1)^3-z(i)^3); %Nm
ABD = [A B; A D];
end
%A_samples(:,:,k) = A;
A_test = A;
%invA_samples(:,:,k) = inv(A);
invA_test = inv(A);
E_x = inv(h*invA_test(1,1));


% Define stress ranges for plotting
sigma_range = linspace(-sigma_max, sigma_max, 101);
N_range = sigma_range*h/100;

% Calculate biaxial stress failure envelopes for Puck and Max Stress criteria
for i = 1:length(sigma_range)
    for j = 1:length(sigma_range)
        sigma1 = sigma_range(i);
        sigma2 = sigma_range(j);
        maxstress_envelope1(i,j) = 0;
        puck_envelope1(i,j) = 0; 
        % Maximum Stress citerion
        if sigma1 <=0 && abs(sigma1/X_C) > 1 % in compression or tension sigma always positive
            maxstress_envelope1(i,j) = 1;
        elseif sigma1 > 0 && abs(sigma1/X_T) > 1
            maxstress_envelope1(i,j) = 1;
        elseif sigma2 <=0 && abs(sigma2/Y_C) > 1 % in compression or tension sigma always positive
            maxstress_envelope1(i,j) = 1;
        elseif sigma2 > 0 && abs(sigma2/Y_T) > 1
            maxstress_envelope1(i,j) = 1;
        end
        %Biaxial Maximum Stress (Puck - Fiber Failure)
        fe = fiberfailure(sigma1,sigma2,0,'c',X_T,X_C,E_x,E1_f,nu12);
        if fe >= 1
            puck_envelope1(i,j) = 1;
        end
    end
end


% Find FPF and LPF for Puck and Max Stress criteria
for i = 1:length(N_range)
    for j = 1:length(N_range)
        F = [N_range(i);N_range(j);0]; % looped to creat different biaxial forces
        strain_glo = invA_test*F;
        max_fe(i,j) = 0;
        stress_glo = Qbar*strain_glo; % global sigmaxx

        for l = 1:Nplies
            % Calculations of Strains and Stresses
            [eps_loc] = strain_gtol(strain_glo,thetadb(l));
            [sigma_loc] = stress_gtol(stress_glo,thetadb(l));% ply i angle in radians, from bottom
            sigma_loc = Q*eps_loc;
            
            % Failure index with Maximum Stress Criterion
            [FI_1(i,j,l),FI_2(i,j,l),FI_3(i,j,l)]= MaxStress(sigma_loc,X_T,X_C,Y_T,Y_C,g_12t);
            FI = [FI_1(i,j,l),FI_2(i,j,l),FI_3(i,j,l)];
            max_FI_1 = max(FI_1(:,:,l))
            %
            % isolate FI_1 and FI_2 which are the biggest with respect to
            % each ply so analyse all the FI_1 for each ply find the
            % maximum and divide the sigma1 by this index. do the same for
            % FI_2 and sigma2
            %
            
            % Failure index with Pucks Criterion
            fe(i,j,l) = fiberfailure(sigma_loc(1),sigma_loc(2),sigma_loc(3),'c',X_T,X_C,E_x,E1_f,nu12);  % (num_samples,layer)
            % Calculation of Ply failure
            if fe(i,j,l) > max_fe(i,j)
                max_fe(i,j)= fe(i,j,l); %maximum failure index with relative forces
                numberply(i,j,l) = l;
            end


        end

        

        % norm of force for first ply failure 
        F_norm_puck(i,j) = norm(F./max_fe(i,j));
        F_norm_MaxStress(i,j) = norm(F./max_FI(i,j));

        Force_failure_maxStress(:,i,j) = [F(1)/max_FI(i,j);F(2)/max_FI(i,j);0];
        strain_failure_maxStress(:,i,j) = invA_test*Force_failure_maxStress(:,i,j);

        Force_failure_Puck(:,i,j) = [F(1)/max_FI(i,j);F(2)/max_fe(i,j);0];
        strain_failure_Puck(:,i,j) = invA_test*Force_failure_Puck(:,i,j);
        
        % degradation
        if max_FI(i,j) < 1
            E_temp(i,j) = 0.1*E2;
        end

    end
end

%%


%%





% % Plot biaxial stress failure envelopes for Puck and Max Stress criteria
% figure(1);
% hold on;
% contour(sigma_range, sigma_range, maxstress_envelope1, [0 1], 'g');
% contour(sigma_range,sigma_range, puck_envelope1, [0 1], 'b');
% xlabel('\sigma_1 (MPa)');
% ylabel('\sigma_2 (MPa)');
% title('Biaxial Stress Failure Envelopes');
% legend('Maximum Stress', 'Puck');
% 
% % Plot FPFfor Puck and Max Stress criteria
% figure(2);
% [X, Y] = meshgrid(N_range, N_range);
% 
% % Plot the 3D surface
% surf(X, Y, max_fe);
% 
% % Add labels and title
% xlabel('x');
% ylabel('y');
% zlabel('fe');
% title('3D Surface Plot of fe');
% 
% 
% figure(3);
% % Plot the 3D surface
% surf(X, Y, max_FI);
% 
% % Add labels and title
% xlabel('x');
% ylabel('y');
% zlabel('FI max');
% title('3D Surface Plot of F');
% 
% figure(4);
% % Plot the 3D surface
% surf(X, Y, F_norm_puck);
% 
% % Add labels and title
% xlabel('x');
% ylabel('y');
% zlabel('F_norm');
% title('3D Surface Plot of F Puck');
% 
% figure(5);
% % Plot the 3D surface
% surf(X, Y, F_norm_MaxStress);
% 
% % Add labels and title
% xlabel('x');
% ylabel('y');
% zlabel('F_norm');
% title('3D Surface Plot of F maxStres');




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


