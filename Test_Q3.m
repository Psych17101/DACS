clear; format;
format short g;

% Define material properties and laminate layup
E1 = 135e9;  % Young's modulus in 1-direction (Pa)
E2 = 10.4e9; % Young's modulus in 2-direction (Pa)
G12 = 4.14e9; % Shear modulus (Pa)
sigma_y = 1130e6; % Yield strength (Pa)
sigma_u = 1220e6; % Ultimate strength (Pa)
t = 0.2e-3; % Laminate thickness (m)
theta = [0 90 -45 45]; % Laminate ply angles (degrees)
h = repmat(t/4, 1, 4); % Ply thicknesses (m)

% Define input variables and Monte Carlo parameters
N = 10000; % Number of simulations
Nx_mean = 100; % Mean axial load level (N/mm)
Ny_mean = 100; % Mean transverse load level (N/mm)
Nx_std = 20; % Standard deviation of axial load level
Ny_std = 20; % Standard deviation of transverse load level

% Perform Monte Carlo simulation
Nx = normrnd(Nx_mean, Nx_std, N, 1); % Generate N samples of axial load
Ny = normrnd(Ny_mean, Ny_std, N, 1); % Generate N samples of transverse load
theta_rad = deg2rad(theta); % Convert laminate ply angles to radians
Q = calculateQ(E1, E2, G12, theta_rad); % Calculate the laminate stiffness matrix
[~, FPF] = calculateTsaiWu(sigma_y, sigma_u, Q, h, Nx, Ny); % Calculate Tsai-Wu FPF
[LaminateStrength, ~] = calculateLaminateStrength(sigma_y, sigma_u, Q, h, Nx, Ny); % Calculate laminate strength
P_FPF = sum(FPF)/N; % Calculate the probability of FPF
P_LPF = sum(LaminateStrength)/N; % Calculate the probability of LPF

% Plot results
figure
plot(Nx, Ny, 'b.', 'MarkerSize', 2); % Plot load points
hold on
plot([0, Nx_mean+Nx_std*3], [0, Ny_mean+Ny_std*3], 'r--'); % Plot 3-sigma boundary
xlabel('Nx (N/mm)')
ylabel('Ny (N/mm)')
title('Biaxial Load Space')
axis equal

figure
histogram(FPF, 50, 'Normalization', 'probability'); % Plot FPF histogram
hold on
line([1 1]*sigma_y, ylim, 'Color', 'r', 'LineWidth', 2); % Plot yield strength
line([1 1]*sigma_u, ylim, 'Color', 'g', 'LineWidth', 2); % Plot ultimate strength
xlabel('Tsai-Wu Failure Index')
ylabel('Probability Density')
title(['P(FPF) = ' num2str(P_FPF)])

figure
histogram(LaminateStrength, 50, 'Normalization', 'probability'); % Plot LPF histogram
hold on
