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
sigma_max = 2*10^9; % 500 MPa

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
sigma_range = linspace(-sigma_max, sigma_max, 101); %N\m^2
% Definition of on distributed biaxial load
N_range = sigma_range*h; %N/m

% For each stress index in range increment.
for i = 1:length(sigma_range)
    for j = 1:length(sigma_range)
        % Biaxial stress
        sigma1 = sigma_range(i);
        sigma2 = sigma_range(j);
        F = [N_range(i);N_range(j);0];
        % Calcuation of global strain for first ply failure
        strain_glo = A\F;
        % Calculation of global stresses
        stress_glo = Qbar*strain_glo; % global sigmaxx etc...

        %Initialisation of Failure criterion envelope
        maxstress_envelope2(i,j) = 0;
        puck_envelope2(i,j) = 0; 

        % Calculations of Failure Criterion for each ply 
        for l = 1:Nplies
            % Calculations of local Strains and Stresses
            [eps_loc] = strain_gtol(strain_glo,thetadb(l));
            %[sigma_loc] = stress_gtol(stress_glo,thetadb(l));% ply i angle in radians, from bottom
            sigma_loc = Qbar*eps_loc;
            
            
            % Failure index with Maximum Stress criterion
            [FI_1(i,j,l),FI_2(i,j,l),FI_3(i,j,l)]= MaxStress(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t);
            
            if FI_1(i,j,l) > 1 || FI_2(i,j,l)>1 || FI_3(i,j,l)>1
                maxstress_envelope2(i,j,l) = 1;
            end

            % Failure index with Pucks Criterion
            [FF(i,j,l),IFF(i,j,l)] = PuckCriterion(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t,nu12,E1,'c','n');
            % Calculation of Ply failure
            if FF(i,j,l) > 1 || IFF(i,j,l) >1
                puck_envelope2(i,j,l) = 1;
            end

        end

    end
end

figure;
hold on
contour(N_range, N_range, maxstress_envelope2(:,:,1), [0 1],'r');
contour(N_range, N_range, maxstress_envelope2(:,:,2), [0 1]),'g';
contour(N_range, N_range, maxstress_envelope2(:,:,3), [0 1],'b'); % 45 = -45 Max stress criterion
xlabel('N_x [N]');
ylabel('N_y [N]');
title('Biaxial Stress Failure Envelopes - Maximum Stress');
legend('Maximum Stress 0 ply', 'Maximum Stress 90 ply','Maximum Stress +-45 ply');
figure;
hold on
contour(N_range, N_range,  puck_envelope2(:,:,1), [0 1],'m');
contour(N_range, N_range,  puck_envelope2(:,:,2), [0 1],'y');
contour(N_range, N_range,  puck_envelope2(:,:,3), [0 1],'k');
xlabel('N_x [N]');
ylabel('N_y [N]');
title('Biaxial Stress Failure Envelopes - Puck');
legend('Puck 0 ply', 'Puck 90 ply','Puck +- 45 ply');


%% Find FPF from Puck and Max Stress criteria
%% Increasing for in Nx Direction
fprintf(" Nx Direction Calculation \n")
for i = 1:length(N_range)
        % Initialisation of loading criterion (Nx only)
        F = [N_range(i);0;0]; 
        % Calcuation of global strain for first ply failure
        strain_glo = A\F;
        % Calculation of global stresses
        stress_glo = Qbar*strain_glo; % global sigmaxx etc...
     
        % Calculations of Failure Criterion for each ply 
         for l = 1:Nplies
            % Calculations of local Strains and Stresses
            [eps_loc] = strain_gtol(strain_glo,thetadb(l));
            [sigma_loc] = stress_gtol(stress_glo,thetadb(l));% ply i angle in radians, from bottom
            disp(sigma_loc)
            % Failure index with Maximum Stress criterion
            [FI_FPF1x(i,l),FI_FPF2x(i,l),FI_FPF3x(i,l)]= MaxStress(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t);
            
            %sigma_loc = Qbar*eps_loc; % same as function stress_gtol
            % Failure index with Pucks Criterion
            [FF_FPFx(i,l),IFF_FPFx(i,l)] = PuckCriterion(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t,nu12,E1,'c','n');
         end
       
        % Considering all the plies - find highest index for each pair of 
        max_FI_1(i) = max(FI_FPF1x(i,:));
        max_FI_2(i) = max(FI_FPF2x(i,:));
        max_FF(i) = max(FF_FPFx(i,:));
        max_IFF(i) = max(IFF_FPFx(i,:)); 
       
        % Creation of First Ply Failure Force
        % Maximum Stress Criterion
        if max_FI_1(i) > max_FI_2(i) % Comparing which Failure index is larger
            max_FI(i) = max_FI_1(i);
            ply_failure = find(max_FI_1(i) == FI_FPF1x(i,:));
        else
            max_FI(i) = max_FI_2(i);
            ply_failure = find(max_FI_2(i) == FI_FPF2x(i,:));
        end
        F_FPF_MS1(i) = N_range(i)/max_FI(i);
        
        % Puck Criterion
        if max_FF(i)>max_IFF(i)
            F_FPF_FF(i) = N_range(i)/max_FF(i);
        else 
            F_FPF_FF(i) = N_range(i)/max_IFF(i);
        end
        
        % Puck criterion only considers plys on lamina level - so it

end

fprintf("First Ply Failure for Nx [N] under Maximum Stress criterion Tension = \n")
disp(max(F_FPF_MS1(length(sigma_range))))
fprintf("First Ply Fiber Failure for Nx [N] under Puck Criterion Tension =\n")
disp(max(F_FPF_FF(length(sigma_range))))

fprintf("First Ply Failure for Nx [N] under Maximum Stress criterion Compression = \n")
disp(min(F_FPF_MS1(1)))
fprintf("First Ply Failure for Nx [N] under Puck Criterion Compression =\n")
disp(min(F_FPF_FF(1)))

%% Increasing Force in Ny direction
fprintf("Ny Direction Calculation \n")
for i = 1:length(N_range)
        % Initialisation of loading criterion (Biaxial loading)
        F = [0;N_range(i);0];
        
        % Calcuation of global strain for first ply failure
        strain_glo = A\F;
        max_fe(i) = 0;
        % Calculation of global stresses 
        stress_glo = Qbar*strain_glo; % global sigmaxx etc...

        % Calculations of Failure Criterion for each ply 
         for l = 1:Nplies
            % Calculations of local Strains and Stresses
            [eps_loc] = strain_gtol(strain_glo,thetadb(l));
            [sigma_loc] = stress_gtol(stress_glo,thetadb(l));% ply i angle in radians, from bottom

            % Failure index with Maximum Stress criterion
            [FI_FPF1y(i,l),FI_FPF2y(i,l),FI_FPF3y(i,l)]= MaxStress(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t);

            % Failure index with Pucks Criterion
            [FF_FPFy(i,l),IFF_FPFy(i,l)] = PuckCriterion(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t,nu12,E1,'c','n');
         end
      
        % Considering all the plies - find highest index for each pair of 
        max_FI_1(i) = max(FI_FPF1y(i,:));
        max_FI_2(i) = max(FI_FPF2y(i,:)); % Strength Ratio
        max_FF(i) = max(FF_FPFy(i,:)); % Strength Ratio
        max_IFF(i) = max(IFF_FPFy(i,:)); % Strength Ratio

        % Creation of First Ply Failure Force
        % Maximum Stress Criterion
        if max_FI_1(i) > max_FI_2(i) % Comparing which Failure index is larger
            max_FI(i) = max_FI_1(i);
            ply_failure = find(max_FI_1(i) == FI_FPF1y(i,:));
        else
            max_FI(i) = max_FI_2(i);
            ply_failure = find(max_FI_2(i) == FI_FPF2y(i,:));
        end
        F_FPF_MS2(i) = N_range(i)/max_FI(i);
        
        % Puck Criterion
        if max_FF(i) > max_IFF(i)
            % Fiber fails first 
            F_FPF_FF(i) = N_range(i)/max_FF(i);
        else 
            % Inter fiber fails first
            F_FPF_FF(i) = N_range(i)/max_IFF(i);
        end
end

fprintf("First Ply Failure for Ny [N] under Maximum Stress criterion Tension = \n")
disp(max(F_FPF_MS2(length(sigma_range))))
fprintf("First Ply Fiber Failure for Ny [N] under Puck Criterion Tension =\n")
disp(max(F_FPF_FF(length(sigma_range))))

fprintf("First Ply Failure for Ny [N] under Maximum Stress criterion Compression = \n")
disp(min(F_FPF_MS2(1)))
fprintf("First Ply Failure for Ny [N] under Puck Criterion Compression =\n")
disp(min(F_FPF_FF(1)))


% Something wrong with Ny - What Fiber failure index or interfiber failure
% index do i use for 






%% LPF

% Calculations of Failure index
i = 1;
% Initialisation of 1 iteration of Modulus
E_temp = E2;
iter = 0;
ply_index = 0;
ply_failure = zeros(1,Nplies);
F = [1;0;0];
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
    strain_glo = A/F;
    % Calculation of global stresses
    stress_glo = Qbar(:,:,1)*strain_glo; % global sigmaxx etc.

    for l = 1:Nplies
        % Calculations of local Strains and Stresses
        [eps_loc] = strain_gtol(strain_glo,thetadb(l));
        [sigma_loc] = stress_gtol(stress_glo,thetadb(l));% ply i angle in radians, from bottom

        [FI_LPF_1(i,l),FI_LPF_2(i,l),FI_LPF_3(i,l)]= MaxStress(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t);

        [FF_LPF1(i,l),IFF_LPF1(i,l)] = PuckCriterion(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t,nu12,E1,'c','n');
    end

    if ply_failure == ones(1,Nplies)
        fprintf('LFP1 =')
        LPF1 = F;
        disp(LPF1)
        break
    end

    temp2 = max(FI_LPF_2(~ply_failure));
    temp1 = max(FI_LPF_1(~ply_failure));

    if  temp1 > temp2
        ply_index = find(FI_LPF_1(i,:) == temp1);
    else
        ply_index = find(FI_LPF_2(i,:) == temp2);
        E_temp = 0.1*E_temp;
    end

    iter = iter + 1;
end


% Calculations of Failure index
i = 1;
% Initialisation of 1 iteration of Modulus
E_temp = E2;
iter = 0;
ply_index = 0;
ply_failure = zeros(1,Nplies);
F = [0;1;0];
temp2 = 0;
% 
while E_temp ~= 0 && iter < 50 % While the failure index of our laminate is less than 1
    z = 0:h_ply:h;
    if temp2 < 10^(-25)
        F = F; % looped to creat different biaxial forces
    else
        F = F/temp2;
    end

    FI_1 = zeros(1,Nplies);
    FI_2 = zeros(1,Nplies);
    FI_3 = zeros(1,Nplies);

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
    E_x = inv(h*invA_test(1,1));

    % Calcuation of global strain for first ply failure
    strain_glo = invA_test*F;
    max_fe(i) = 0;
    % Calculation of global stresses
    stress_glo = Qbar(:,:,1)*strain_glo; % global sigmaxx etc.

    for l = 1:Nplies
        % Calculations of local Strains and Stresses
        [eps_loc] = strain_gtol(strain_glo,thetadb(l));
        [sigma_loc] = stress_gtol(stress_glo,thetadb(l));% ply i angle in radians, from bottom
        %sigma_loc = Qbar(:,:,1)*eps_loc;


        %

        % Failure index with Maximum Stress Criterion
        % For each pair of N(i) and N(j) there will be a maximum
        % failure criterion - from all the calculations of each failure
        % index, find the maximum one.

        [FI_1(i,l),FI_2(i,l),FI_3(i,l)]= MaxStress(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t);

        % Failure index with Pucks Criterion
        fe(i,l) = fiberfailure(sigma_loc(1),sigma_loc(2),sigma_loc(3),'c',X_T,X_C,E_x,E1_f,nu12);  % (num_samples,layer)
        % Calculation of Ply failure
        if fe(i,l) > max_fe(i)
            max_fe(i)= fe(i,l); %maximum failure index with relative forces
            numberply(i,l) = l;
        end
    end

    if ply_failure == ones(1,Nplies)
        fprintf('LFP2 =')
        LPF2 = F;
        disp(LPF2)
        break
    end

    temp2 = max(FI_2(~ply_failure));
    temp1 = max(FI_1(~ply_failure));

    if  temp1 > temp2
        ply_index = find(FI_2(i,:) == temp2);
    else
        ply_index = find(FI_1(i,:) == temp1);
        E_temp = 0.1*E_temp;
    end

    iter = iter + 1;
end


%% Plots
% Plot biaxial stress failure envelopes for Puck and Max Stress criteria
figure;
hold on;
div = 10^6;
contour(sigma_range, sigma_range, maxstress_envelope1, [0 1], 'g');
contour(sigma_range,sigma_range, puck_envelope1, [0 1], 'b');
scatter(F_FPF_MS1(1)/h,0,'x') % Nx FPF in Compression Max stress
scatter(F_FPF_MS1(1000)/h,0,'x') % Nx FPF in Tension Max stress
scatter(F_FPF_FF1(1)/h,0,'o') % Nx FPF in Compression Puck
scatter(F_FPF_FF1(1000)/h,0,'o') % Nx FPF in Tension Puck
scatter(LPF1(1)/h,0,'*') % Last ply Failure in Tension
% LFP in compression need
scatter(0,F_FPF_MS2(1)/h,'x') % Ny FPF in Compression Max stress
scatter(0,F_FPF_MS2(1000)/h,'x') % Ny FPF in Tension Max Stress
scatter(0,F_FPF_FF2(1)/h,'o') % Ny FPF in Compression Puck
scatter(0,F_FPF_FF2(1000)/h,'o') % Ny FPF in Tension Puck
scatter(0,LPF2(2)/h,'*')
xlabel('\sigma_1 (Pa)');
ylabel('\sigma_2 (Pa)');
title('Biaxial Stress Failure Envelopes');
legend('Maximum Stress', 'Puck', ...
    'FPF Maximum stress -Nx','FPF Maximum stress +Nx', ...
    'FPF Puck Criterion -Nx','FPF Puck Criterion +Nx', ...
    'LPF Maximum stress Criterion Nx', ...
    'FPF Maximum stress -Ny','FPF Maximum stress +Ny', ...
    'FPF Puck Criterion -Ny','FPF Puck Criterion +Ny', ...
    'LPF Maximum stress Criterion Ny');


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

function [FF,IFF] = PuckCriterion(sigma1,sigma2,sigma3,X_T,X_C,Y_T,Y_C,G12_t,nu12,E1,fiber,print)
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
    if print == 'y'
        fprintf("fiber failure Tension\n")
    end
else 
    FF = 1/(-X_C)*(sigma1 - (nu12-nu12*m*E1/Ex)*(sigma2+sigma3));
    if print == 'y'
        fprintf("Fibre failure Compression\n")
    end
end

if sigma2 > 0  && Y_T > sigma2 %Condition Mode A
    R_Tll = G12_t;
    R_Tt = Y_T;
    A_ = (1/R_Tt - ptTll/R_Tll)*sigma2;
    B = sigma3/R_Tll;
    C = ptTll/R_Tll*sigma2;
    IFF = sqrt(A_^2 + B^2) + C;
    if print == 'y'
        fprintf("Interfiber failure Mode A\n")
    end
elseif sigma2 < 0 && -R_TTA < sigma2  %  Condition Mode B
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
stress_loc=R*T*inv(R)*stress_glo;
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


