clear; format;
close;
format short g;

% Define material properties
E1   = 165.22*10^9 ;% Pa - direction modulus of the lamina
E2   = 8.4443*10^9  ; % Pa
nu12 = .35218 ;
nu21 = 0.0185;
G12  = 6.444*10^9  ; % Pa
X_T = 1923.7*10^6;
X_C = 1480*10^6;
Y_T = 107*10^6;
Y_C = 220*10^6;
g_12t = 152.4*10^6;

%% Initialisation of Composites Model
% Laminate definition (plies of equal thickness)

thetadt = [0 90 +45 -45 -45 +45 90 0 0 90 +45 -45 -45 +45 90 0];% ply angles in degrees, from top??? Not bottom?
Nplies  = length(thetadt);
thetadb = fliplr(thetadt); % ply angles in degrees, from bottom
h_ply   = 127.84*10^(-6);   % SI units, meters
h       = Nplies * h_ply ;

z       = -h/2:h_ply:h/2;

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
for i = 1:Nplies
    [Qbar(:,:,i),Sbar(:,:,i)] = QbarandSbar(thetadb(i),moduli);
    A = A + Qbar(:,:,i) * (z(i+1)-z(i)) ; %N/m, 
    B = B + (1/2)* Qbar(:,:,i)* (z(i+1)^2-z(i)^2); %N
    D = D + (1/3)* Qbar(:,:,i) * (z(i+1)^3-z(i)^3); %Nm
    ABD = [A B; B D];
end
% Initialisation of Equivalent Modulus of Laminate
invA = inv(A);
E_x = inv(h*invA(1,1));

%% Plot and Calculations biaxial stress failure envelopes for Puck and Max Stress criteria
% Define stress ranges for plotting
a = 0.01;
sigma_range = linspace(-a*Y_C, a*Y_T, 31); %N\m^2

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
        stress_glo = Qbar(:,:,1)*strain_glo; % global sigmaxx etc...

        % Calculations of Failure Criterion for each ply
        for l = 1:Nplies
            % Calculations of local Strains and Stresses
            [eps_loc] = strain_gtol(strain_glo,thetadb(l));
            [sigma_loc] = stress_gtol(stress_glo,thetadb(l));% ply i angle in radians, from bottom

            % Failure index with Maximum Stress criterion
            [FI_1(i,j,l),FI_2(i,j,l),FI_3(i,j,l)]= MaxStress(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t);

            % Failure index with Pucks Criterion
            [FF(i,j,l),IFF(i,j,l)] = PuckCriterion(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t,nu12,E1,'c','n');

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
        FPF_MS(i,j,:) = F/max_FI_MS(i,j);
        Nx_FPF_MS(i,j) = FPF_MS(i,j,1);
        Ny_FPF_MS(i,j) = FPF_MS(i,j,2);

        % Puck Criterion
        FPF_PUCK(i,j,:) = F/max_FI_PUCK(i,j);
        Nx_FPF_PUCK(i,j) = FPF_PUCK(i,j,1);
        Ny_FPF_PUCK(i,j) = FPF_PUCK(i,j,2);
    end
end


%% LPF - MS
for i = 1:length(sigma_range)
    for j = 1:length(sigma_range)
        % Calculations of Failure index
        % Initialisation of 1 iteration of Modulus
        iter = 0;
        ply_index = 0;
        ply_failure = zeros(1,Nplies);
        F_LPF = [N_range(i);N_range(j);0];
        %disp(F_LPF)
        max_FI_MS_LPF(i,j) = 0;

        while iter < 50 
            z = -h/2:h_ply:h/2;
            if max_FI_MS_LPF(i,j) == 0
                F_LPF(:) = F_LPF(:); % looped to creat different biaxial forces
            else
                F_LPF(:) = F_LPF(:)/max_FI_MS_LPF(i,j);
            end

            FI_LPF_1(i,j,:) = zeros(1,Nplies);
            FI_LPF_2(i,j,:) = zeros(1,Nplies);
            FI_LPF_3(i,j,:) = zeros(1,Nplies);

            % ABD reset
            A = zeros(3,3);
            B = zeros(3,3);
            D = zeros(3,3);
            Qbar = zeros(3,3);
            Sbar = zeros(3,3);

            % Calculation of Stiffness Matrix
            [S, Q] = ReducedComplianceStiffness(E1,E2,nu12,G12);
            [moduli]= [E1 E2 nu12 G12];

            zbar = zeros(1, Nplies);
            for l = 1:Nplies
                zbar(l) = -(h + h_ply)/2 + l*h_ply;
                % For each ply we calculate the ABD Matrix
                [Qbar(:,:,l),Sbar(:,:,l)] = QbarandSbar(thetadb(l),moduli);

                if ply_failure(l) == 2
                    Qbar(2,2,l) = 0.1*Qbar(2,2,l); % Degrading properties if matrix failure
                elseif ply_failure(l) == 1
                    Qbar(:,:,l) = zeros(3,3); % Zero properties within the laminate
                end

                A = A + Qbar(:,:,l) * (z(l+1)-z(l)) ; %N/m, right dimensions?
                B = B + (1/2) * Qbar(:,:,l) * (z(l+1)^2-z(l)^2); %N
                D = D + (1/3) * Qbar(:,:,l) * (z(l+1)^3-z(l)^3); %Nm
                ABD = [A B; B D];
            end

            % Calcuation of global strain for failure
            strain_glo = inv(A)*F_LPF;
           
            % Calculation of global stresses
            stress_glo = Qbar(:,:,1)*strain_glo; % global sigmaxx etc.

            for l = 1:Nplies
                % Calculations of local Strains and Stresses
                [eps_loc] = strain_gtol(strain_glo,thetadb(l));
                [sigma_loc] = stress_gtol(stress_glo,thetadb(l));% ply i angle in radians, from bottom
                [FI_LPF_1(i,j,l),FI_LPF_2(i,j,l),FI_LPF_3(i,j,l)]= MaxStress(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t);
            end

            if ply_failure ~= zeros(1,Nplies) % Checking if all plies fail
                Nx_LPF_MS(i,j) = F_LPF(1);
                Ny_LPF_MS(i,j) = F_LPF(2);
                break
            end

            % Considering all the plies - find highest failure index of loading ration (i,j)
            % Maximum Stress
            max_FI_1_LPF(i,j) = max(FI_LPF_1(i,j,~ply_failure)); % Maximum FI_1 over all plies for loading ratio (i,j)
            max_FI_2_LPF(i,j) = max(FI_LPF_2(i,j,~ply_failure));
            max_FI_3_LPF(i,j) = max(FI_LPF_3(i,j,~ply_failure));
            % Highest Failure index according to maxmimum Stress
            max_FI_MS_LPF(i,j) = max([max_FI_1_LPF(i,j),max_FI_2_LPF(i,j),max_FI_3_LPF(i,j)]);

            if max_FI_MS_LPF(i,j) == max_FI_1_LPF(i,j) % if fibre failure
                ply_index = find(max_FI_MS_LPF(i,j) == FI_LPF_1(i,j));
                for n = 1:length(ply_index) % Detection of ply failure
                    if ply_index ~= 0
                        ply_failure(ply_index(n)) = 1;
                    end
                end
            elseif max_FI_MS_LPF(i,j) == max_FI_2_LPF(i,j) % if Matrix failure
                ply_index = find(FI_LPF_2(i,j,:) == max_FI_MS_LPF(i,j,:));
                for n = 1:length(ply_index) 
                    if ply_index ~= 0
                        ply_failure(ply_index(n)) = 2;
                    end
                end
            else
                ply_index = find(max_FI_MS_LPF(i,j) == FI_LPF_3(i,j,:)); 
                for n = 1:length(ply_index)
                    if ply_index ~= 0
                        ply_failure(ply_index(n)) = 2;
                    end
                end
            end

            iter = iter + 1;
        end
    end
end

%% LPF - PUCK
% Calculations of Failure index
        % Initialisation of 1 iteration of Modulus
for i = 1:length(sigma_range)
    for j = 1:length(sigma_range)
        iter = 1;
        ply_index = 0;
        ply_failure = zeros(1,Nplies);
        F_LPF = [N_range(i);N_range(j);0];
        max_FI_PUCK_LPF(i,j) = 0;

        while iter < 50 % While the failure index of our laminate is less than 1
            z = -h/2:h_ply:h/2;
            if max_FI_PUCK_LPF(i,j) == 0
                F_LPF(:) = F_LPF(:); % looped to creat different biaxial forces
            else
                F_LPF(:) = F_LPF(:)/max_FI_PUCK_LPF(i,j);
            end
            FF_LPF(i,j,:) = zeros(1,Nplies);
            IFF_LPF(i,j,:) = zeros(1,Nplies);

            % ABD reset
            A = zeros(3,3);
            B = zeros(3,3);
            D = zeros(3,3);
            Qbar = zeros(3,3);
            Sbar = zeros(3,3);

            % Calculation of Stiffness Matrix
            [S, Q] = ReducedComplianceStiffness(E1,E2,nu12,G12);
            [moduli]= [E1 E2 nu12 G12];

            zbar = zeros(1, Nplies);

            for l = 1:Nplies
                zbar(l) = -(h + h_ply)/2 + l*h_ply;
                % For each ply we calculate the ABD Matrix
                [Qbar(:,:,l),Sbar(:,:,l)] = QbarandSbar(thetadb(l),moduli);

                if ply_failure(l) == 2
                    Qbar(2,2,l) = 0.1*Qbar(2,2,l);
                elseif ply_failure(l) == 1
                    Qbar(:,:,l) = zeros(3,3);
                end

                A = A + Qbar(:,:,l) * (z(l+1)-z(l)) ; %N/m
                B = B + (1/2) * Qbar(:,:,l) * (z(l+1)^2-z(l)^2); %N
                D = D + (1/3) * Qbar(:,:,l) * (z(l+1)^3-z(l)^3); %Nm
                ABD = [A B; B D];
            end

            % Calcuation of global strain for first ply failure
            strain_glo = inv(A)*F_LPF;
            % Calculation of global stresses
            stress_glo = Qbar(:,:,1)*strain_glo; % global sigmaxx etc.

            for l = 1:Nplies
                % Calculations of local Strains and Stresses
                [eps_loc] = strain_gtol(strain_glo,thetadb(l));
                [sigma_loc] = stress_gtol(stress_glo,thetadb(l));% ply i angle in radians, from bottom

                if ply_failure(l) == 1
                    [FF_LPF(i,j,l),IFF_LPF(i,j,l)]= PuckCriterion(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t,nu12,0,'c','n');
                else
                    [FF_LPF(i,j,l),IFF_LPF(i,j,l)]= PuckCriterion(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,g_12t,nu12,E1,'c','n');
                end
                
            end

            if ply_failure ~= zeros(1,Nplies)
                Nx_LPF_PUCK(i,j) = F_LPF(1);
                Ny_LPF_PUCK(i,j) = F_LPF(2);
                break
            end

            % Considering all the plies - find highest failure index of loading ration (i,j)
            % PUCK
            max_FF_LPF(i,j) = max(FF_LPF(i,j,~ply_failure)); % Maximum FI_1 over all plies for loading ratio (i,j)
            max_IFF_LPF(i,j) = max(IFF_LPF(i,j,~ply_failure));
            % Highest Failure index according to maxmimum Stress
            max_FI_PUCK_LPF(i,j) = max([max_FF_LPF(i,j),max_IFF_LPF(i,j)]);

            if max_FI_PUCK_LPF(i,j) == max_FF_LPF(i,j) % if fibre failure
                ply_index = find(max_FI_PUCK_LPF(i,j) == FF_LPF(i,j,:));
                for n = 1:length(ply_index)
                    if ply_index ~= 0
                        ply_failure(ply_index(n)) = 1;
                    end
                end
            elseif max_FI_PUCK_LPF(i,j) == max_IFF_LPF(i,j)
                ply_index = find(max_FI_PUCK_LPF(i,j) == IFF_LPF(i,j,:));
                for n = 1:length(ply_index)
                    if ply_index ~= 0
                        ply_failure(ply_index(n)) = 2;
                    end
                end
            end
  
            iter = iter + 1;
        end
    end
end

%% Plotting

% MS FPF && LPF
vec_LR = LR(:);
[sorted_vecLR sorting_index] = sort(vec_LR,'descend');

Nx_FMS = Nx_FPF_MS(sorting_index);
Ny_FMS = Ny_FPF_MS(sorting_index);
Nx_LMS = Nx_LPF_MS(sorting_index);
Ny_LMS = Ny_LPF_MS(sorting_index);

vec_Nx_FMS = Nx_FMS(:);
vec_Ny_FMS = Ny_FMS(:);

vec_Nx_LMS1 = Nx_LMS(Ny_LMS(:)>-8*10^5);
vec_Ny_LMS1 = Ny_LMS(Ny_LMS(:)>-8*10^5);

vec_Nx_LMS = vec_Nx_LMS1(vec_Ny_LMS1(:)<4*10^5);
vec_Ny_LMS = vec_Ny_LMS1(vec_Ny_LMS1(:)<4*10^5);


P_FMS = [vec_Nx_FMS,vec_Ny_FMS];
k_FMS = convhull(P_FMS);
P_LMS = [vec_Nx_LMS,vec_Ny_LMS];
k_LMS = convhull(P_LMS);

%%
vec_LR = LR(:);
[sorted_vecLR sorting_index] = sort(vec_LR,'descend');
Nx_FP = Nx_FPF_PUCK(sorting_index);
Ny_FP = Ny_FPF_PUCK(sorting_index);
Nx_LP = Nx_LPF_PUCK(sorting_index);
Ny_LP = Ny_LPF_PUCK(sorting_index);

vec_Nx_FP = Nx_FP(:);
vec_Ny_FP = Ny_FP(:);

vec_Nx_LP1 = Nx_LP(Ny_LP(:)>-0.5*10^6);
vec_Ny_LP1 = Ny_LP(Ny_LP(:)>-0.5*10^6);
% vec_Nx_LP1 = Nx_LP(:);
% vec_Ny_LP1 = Ny_LP(:);

vec_Nx_LP = vec_Nx_LP1(vec_Ny_LP1(:)<4*10^5);
vec_Ny_LP = vec_Ny_LP1(vec_Ny_LP1(:)<4*10^5);
% vec_Nx_LP = vec_Nx_LP1(:);
% vec_Ny_LP = vec_Ny_LP1(:);

P_FP = [vec_Nx_FP,vec_Ny_FP];
k_FP = convhull(P_FP);
P_LP = [vec_Nx_LP,vec_Ny_LP];
k_LP = convhull(P_LP);


figure(1);
hold on
grid on
plot(P_FP(:,1),P_FP(:,2),'b*')
plot(P_FP(k_FP,1),P_FP(k_FP,2),'b')
plot(P_LP(:,1),P_LP(:,2),'r*')
plot(P_LP(k_LP,1),P_LP(k_LP,2),'r')
xlabel('$Nx$[N]','interpreter','Latex')
ylabel('$Ny$[N]','interpreter','Latex')
legend('FPF - Puck','FPF envelope - Puck','LPF - Puck','LPF envelope - Puck')
title('Comparison PUCK - FPF & LPF')
exportgraphics(gcf,'PUCK_Comp_points.png','Resolution',500);

% Comparison MS and PUCK
figure(2);
hold on
grid on
plot(P_FP(:,1),P_FP(:,2),'b*')
plot(P_FP(k_FP,1),P_FP(k_FP,2),'b')
plot(P_FMS(:,1),P_FMS(:,2),'r*')
plot(P_FMS(k_FMS,1),P_FMS(k_FMS,2),'r')
xlabel('$Nx$[N]','interpreter','Latex')
ylabel('$Ny$[N]','interpreter','Latex')
legend('FPF - Puck','FPF envelope - Puck','FPF - MS','LPF envelope - MS')
title('Comparison FPF - Puck & MS')
exportgraphics(gcf,'FPF_Comp_points.png','Resolution',500);


%% Envelopes 
% Comparison MS and PUCK FPF

[xcw1, ycw1] = poly2cw(vec_Nx_FMS,vec_Ny_FMS);
K = convhull(xcw1,ycw1);
xy = [xcw1,ycw1];
xy_FMS = xy(K,:);
[xcw2, ycw2] = poly2cw(vec_Nx_FP,vec_Ny_FP);
K = convhull(xcw2,ycw2);
xy = [xcw2,ycw2];
xy_FP = xy(K,:);

N_range = linspace(-6*10^6,1*10^6);
Nx = cosd(30)*400*10^3;
Ny = sind(30)*400*10^3;
Nx_1 = cosd(30)*100*10^3;
Ny_1 = sind(30)*100*10^3;
Ny_range = sind(30)/cosd(30)*N_range;

% Define angle range to remove points
angle_min = -135; % minimum angle
angle_max = -90; % maximum angle
[xcw1, ycw1] = poly2cw(vec_Nx_LMS,vec_Ny_LMS);
K = convhull(xcw1,ycw1);
xy = [xcw1,ycw1];
xy_LMS = xy(K,:);
% Remove identified points from polygon xy1
angles = atan2d(xy_LMS(:,2), xy_LMS(:,1));
angles = mod(angles + 180, 360) - 180;
idx = (angles >= angle_min) & (angles <= angle_max);
xy_LMS(idx, :) = [];
P_LMS(idx,:) = [];

angle_min = 45; % minimum angle
angle_max = 90; % maximum angle
[xcw2, ycw2] = poly2cw(vec_Nx_LP,vec_Ny_LP);
xy_LP = [xcw2,ycw2];

% Remove identified points from polygon xy2
angles = atan2d(xy_LP(:,2), xy_LP(:,1));
angles = mod(angles + 180, 360) - 180;
idx = (angles >= angle_min) & (angles <= angle_max);
xy_LP(idx, :) = [];
K = convhull(xy_LP(:,1),xy_LP(:,2));
xy_LP = xy_LP(K,:);
P_LP(idx,:) = [];

figure(3); % PLY FAILURE ITERATION
hold on
grid on
axis equal
scatter(Nx,Ny,'go')
scatter(Nx_1,Ny_1,'ro')
plot(P_FMS(:,1),P_FMS(:,2),'r*')
plot(xy_FMS(:,1), xy_FMS(:,2), 'g')
plot(P_FP(:,1),P_FP(:,2),'b*')
plot(xy_FP(:,1), xy_FP(:,2), 'k')
plot(N_range,Ny_range,'b-')
xlabel('$Nx$[N]','interpreter','Latex')
ylabel('$Ny$[N]','interpreter','Latex')
legend('Loading: 400 N/m','Loading: 100 N/m','FPF envelope - Puck','FPF envelope - MS','Slope')
title('Comparison FPF - Puck & MS')
exportgraphics(gcf,'FPF_Comp.png','Resolution',500);

figure(4); % Comparison Last ply Failure
hold on
grid on
axis equal
plot(P_LMS(:,1), P_LMS(:,2), 'g*')
plot(P_LP(:,1), P_LP(:,2), 'k*')
plot(xy_LMS(:,1), xy_LMS(:,2), 'g')
plot(xy_LP(:,1), xy_LP(:,2), 'k')
xlabel('$Nx$[N]','interpreter','Latex')
ylabel('$Ny$[N]','interpreter','Latex')
legend('LPF envelope - MS','LPF envelope - PUCK')
title('Comparison LPF - Puck & MS')
exportgraphics(gcf,'LPF_Comp.png','Resolution',500);

figure(5); % Maximum stress criterion
hold on
grid on
axis equal
plot(xy_FMS(:,1), xy_FMS(:,2), 'b')
plot(xy_LMS(:,1), xy_LMS(:,2), 'r')
xlabel('$Nx$[N]','interpreter','Latex')
ylabel('$Ny$[N]','interpreter','Latex')
legend('FPF envelope - MS','LPF envelope - MS')
title('Comparison Maximum Stress - FPF & LPF')
exportgraphics(gcf,'MS_Comp.png','Resolution',500);

figure(6); % Pucks Criterion
hold on
grid on
axis equal
plot(xy_FP(:,1), xy_FP(:,2), 'b')
plot(xy_LP(:,1), xy_LP(:,2), 'r')
xlabel('$Nx$[N]','interpreter','Latex')
ylabel('$Ny$[N]','interpreter','Latex')
legend('FPF envelope - PUCK','LPF envelope - PUCK')
title('Comparison Puck - FPF & LPF')
exportgraphics(gcf,'PUCK_Comp.png','Resolution',500);

figure(7); % Pucks Criterion
hold on
grid on
axis equal
plot(xy_FMS(:,1), xy_FMS(:,2), 'b')
plot(xy_FP(:,1), xy_FP(:,2), 'r')
xlabel('$Nx$[N]','interpreter','Latex')
ylabel('$Ny$[N]','interpreter','Latex')
legend('FPF envelope - PUCK','LPF envelope - PUCK')
title('Comparison FPF')
exportgraphics(gcf,'FPF_Comp2.png','Resolution',500);


%% Functions
function [FF,IFF] = PuckCriterion(sigma1,sigma2,sigma3,X_T,X_C,Y_T,Y_C,G12_t,nu12,E1,fiber,print)
% This function returns the Failure indices under the PUCK criterion
% for fiber-reinforced materials.
% There are inputs are the local stresses with the materials properties and strengths
% the fiber input and returns each fiber failure and inter fibre failure
% index

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
    Ex = 270e9; % GPa  
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
% This function returns the Failure indices under the maximum stress criterion
% for fiber-reinforced materials.
% There are inputs are the local stresses with the materials strengths
% and returns each failure index
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
% This function returns the transformed reduced stiffness
% matrix for fiber-reinforced materials.
% There are four arguments representing four
% material constants. The size of the reduced
% stiffness and compliance matrix is 3 x 3.
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

