clear; format;
close all;
format short g;

%% Notes
% I am a bit stuck on the creation of the probability
% Currently taking random values of material propoerties
% Not sure if calculation of first ply failure is correct (should be) (all values tend to be lower that 1 = no ply failure)
% Not sure about units of force and material properties
% Problems with plotting


%% Initialisationation
layup = [0 +45 -45 90 90 -45 +45 0 ];
Nplies = length(layup);
thetadb = fliplr(layup); % ply angles in degrees, from bottom;
h_ply  = 0.125*10^(-3);           % SI units, meters
h      = Nplies * h_ply ;

z = 0:h_ply:h;

for i = 1:Nplies
  zbar(i) = - (h + h_ply)/2 + i*h_ply;
end
% material properties (mean and standard deviation)
mean_E1 = 140e6; std_E1 = 20e6;
mean_E2 = 10e6; std_E2 = 0.8e6;
mean_G12 = 5e6; std_G12 = 0.3e6;
mean_v12 = 0.2; std_v12 = 0.01;
X_T = 1500.e6; % Fiber
X_C = 1200.e6; %
Y_T = 50.e6; % Matrix
Y_C = 250.e6;
S_t = 70.e6;


% load and orientation
N = [100000, 400000]; % resultant force in N/mm
theta = [30, 30]; % angle of resultant force with respect to X-axis in radians

% Monte Carlo simulation parameters
num_simulations = 100; % number of Monte Carlo simulations to run
num_samples = 30; % number of random samples to generate for each simulation

for i = 1:length(N)
    F = [N(i)*cos(theta(i)); N(i)*sin(theta(i)); 0]; % load vector
    ply_failure = zeros(num_simulations, 1); % array to store ply failure results
    for j = 1:num_simulations
        % generate random samples for material properties
        E1_samples = normrnd(mean_E1, std_E1, 1, num_samples);
        E2_samples = normrnd(mean_E2, std_E2, 1, num_samples);
        G12_samples = normrnd(mean_G12, std_G12, 1, num_samples);
        v12_samples = normrnd(mean_v12, std_v12, 1, num_samples);
        
        % calculate the transformed stiffness matrix for each sample
        
        strain_samples = zeros(3, num_samples);
        stress_samples = zeros(3, num_samples);

        for k = 1:num_samples
            A = zeros(3,3);
            B = zeros(3,3);
            D = zeros(3,3);
            [S, Q] = ReducedComplianceStiffness(E1_samples(k),E2_samples(k),v12_samples(k),G12_samples(k));
            [moduli] = [E1_samples(k),E2_samples(k),v12_samples(k),G12_samples(k)]; 

            for l = 1:Nplies
                % For each ply we calculate the ABD Matrix
                [Qbar,Sbar] = QbarandSbar(thetadb(l),moduli);

                A = A + Qbar * (z(i+1)-z(i)) ; %N/m, right dimensions?
                B = B + (1/2)*Qbar * (z(i+1)^2-z(i)^2); %N 
                D = D + (1/3)*Qbar * (z(i+1)^3-z(i)^3); %Nm
                ABD = [A B; A D];
            end
            A_samples(:,:,k) = A;
            A_test = A;
            invA_samples(:,:,k) = inv(A);
            invA_test = inv(A);
            E_x(k) = inv(h*invA_test(1,1));
            strain_samples(:,k) = invA_samples(:,:,k)*F;
            max_fe(k) = 0;

            for l = 1:Nplies
                % Calculations of Strains and Stresses
                stress_samples(:,k) = Q*strain_samples(:,k); % global sigmaxx
                [eps_loc] = strain_gtol(strain_samples(:,k),thetadb(l)); % ply i angle in radians, from bottom
                sigma_loc = Q*eps_loc;
                fe(k,l) = puckfailure(sigma_loc,E1_samples(k),E_x(k),v12_samples(k),X_T,X_C);  % (num_samples,layer)
                % Calculation of Ply failure 
                if fe(k,l) > max_fe(k)
                    max_fe(k) = fe(k,l);
                    numberply(k,l) = l;
                end
            end    
      
             
        end
        
        % calculate the probability of failure for each ply
        % check if any ply has failed for each sample
        %ply_failure(j) = max_fe > 1;
    end
    
    % calculate the probability of failure for each load and orientation
    %P_f_convergence(k) = sum(ply_failure) / num_simulations;
end

% plot the probability of failure as a function of load and orientation
% figure(1)
% plot(N, P_f)
% title('Probability of First Ply Failure vs. Load and Orientation')
% xlabel('Load (N/mm)')
% ylabel('Probability of Failure')
% legend('0°', '30°', '90°', '-30°', '-90°')

num_simulations_vec = [1,5,10];
% check convergence using the theorem of large numbers
figure(2)
%plot(num_simulations_vec, P_f_convergence)
title('Convergence of Probability of First Ply Failure')
xlabel('Number of Simulations')
ylabel('Probability of Failure')

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

function fe = puckfailure(sigma_loc,E1,Ex,nu12,X_T,X_C)
m = 1.3;
    if sigma_loc(1) > 0
        fe = 1/X_T*(sigma_loc(1) - (nu12-nu12*m*Ex/E1)*(sigma_loc(2)+sigma_loc(3)));
    else
        fe = 1/(-X_C)*(sigma_loc(1) - (nu12-nu12*m*Ex/E1)*(sigma_loc(2)+sigma_loc(3)));
    end
end
