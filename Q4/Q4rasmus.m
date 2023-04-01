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
h_ply  = 0.125*10^(-3);           % SI units, mmm
h      = Nplies * h_ply ;

z = 0:h_ply:h;

for i = 1:Nplies
  zbar(i) = - (h + h_ply)/2 + i*h_ply; % mm
end

% material properties (mean and standard deviation) in Pa
mean_E1 = 165e9; std_E1 = 17e9;
mean_E2 = 8.44e9; std_E2 = 7.5e9;
mean_G12 = 7e9; std_G12 = 0.1e9;
mean_v12 = 0.35; std_v12 = 0.18;
mean_X_T = 1920.e6; std_X_T=108.65*1e8;% Fiber
X_C = 1200.e6; %
mean_Y_T = 107.e6; std_Y_T=9.35*1e6;% Matrix
Y_C = 250.e6;
S_t = 70.e6;


% load and orientation
N = 400*1e3; % resultant force in N/m
% Confused about units
theta = 30; % angle of resultant force with respect to X-axis in radians

% Monte Carlo simulation parameters
num_simulations = 10; % number of Monte Carlo simulations to run
num_samples = 1; % number of random samples to generate for each simulation

F = [N*cosd(theta); N*sind(theta); 0]; % load vector N/m
ply_failure = zeros(num_simulations, 1); % array to store ply failure results
error=1;
pflist=[];
hh = animatedline('Color','r');
axis([0,10000,0,5])
FN=1;
FNlast=1;
while error>1e-3 %run until converged
    % generate random samples for material properties
    
    n=0;
    Truth=1;
    while Truth==1
        E1_samples = normrnd(mean_E1, std_E1, 1); % Pa
        E2_samples = normrnd(mean_E2, std_E2, 1); % Pa
        G12_samples = normrnd(mean_G12, std_G12, 1); % Pa
        v12_samples = normrnd(mean_v12, std_v12, 1);
        X_T = normrnd(mean_X_T, std_X_T, 1);
        Y_T = normrnd(mean_Y_T, std_Y_T, 1);
        % calculate the transformed stiffness matrix for each sample

        %strain_samples = zeros(3);
        %stress_samples = zeros(3);
   
        A = zeros(3,3);
        B = zeros(3,3);
        D = zeros(3,3);
        [S, Q] = ReducedComplianceStiffness(E1_samples,E2_samples,v12_samples,G12_samples);
        [moduli] = [E1_samples,E2_samples,v12_samples,G12_samples]; 

        for l = 1:Nplies
            % For each ply we calculate the ABD Matrix
            [Qbar,Sbar] = QbarandSbar(thetadb(l),moduli);

            A = A + Qbar * (z(i+1)-z(i)) ; %N/mm
            B = B + (1/2)*Qbar * (z(i+1)^2-z(i)^2); % N 
            D = D + (1/3)*Qbar * (z(i+1)^3-z(i)^3); %N mm
            ABD = [A B; B D];
        end


        A_test = A;
        invA_samples = inv(A);
        E_x = inv(h*invA_samples);
        strain_samples = A\F;%invA_samples*F; ;
        max_fe= 0;

        for l = 1:Nplies
            % Calculations of Strains and Stresses
            stress_samples = Q*strain_samples; % global sigmaxx % Pa
            [eps_loc] = strain_gtol(strain_samples(:),thetadb(l)); % ply i angle in radians, from bottom
            sigma_loc = Q*eps_loc;


            % Calculations of Puck criterion for each simulation run,
            % each random variable chosen and each ply 
            fe = puckfailure(sigma_loc,E1_samples,E_x,v12_samples,X_T,X_C);  % (num_samples,layer)

            % Calculation of Ply failure 
            % Calculation of max failure index under puck criterion
                % for each simulation and each random sample
            %if fe > max_fe
            %    max_fe = fe; % failure index under puck criterion
            %    %numberply(j,k,l) = l;
            %end
            if fe>1
                Truth=0;
            end
        
        end  
    n=n+1;
    end
    pflist=[pflist, 1/n];%P_f^R
    % is the failure index bigger than 1?
    FNlast2=FNlast;
    FNlast=FN;
    FN=sum(pflist)/length(pflist);
    
    %plot(length(FN), FN); hold on;
    ROC=abs(FNlast2-FNlast)/abs(FNlast-FN);
    %addpoints(hh,length(pflist),FN);
    
    addpoints(hh, length(pflist), ROC); 
    drawnow;
    error=FN-FNlast;
    
    
end

% calculate the probability of failure for each load and orientation
%P_f_convergence(k) = sum(ply_failure) / num_simulations;

%%
% plot of maximum failure index against number of simulations run
figure(1)
array_num_simulations = 1:num_simulations;
hold on
for i = 1:num_samples
    plot(array_num_simulations,max_fe(:,i))
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
