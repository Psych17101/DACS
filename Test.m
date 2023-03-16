clear; format;
format short g;

% Define material properties
E1   = 140.e9 ; % Pa
nu12 = .3 ;
E2   = 10.e9  ; % Pa
G12  = 5.e9  ; % Pa
nu21 = nu12 * E2 / E1 ;
X_T = 1500.e6;
X_C = 1200.e6;
Y_T = 50.e6;
Y_C = 250.e6;
g_12t = 70.e6;
sigma_max = 90.e8;

% Define failure criteria
tau_12f = g_12t;
puck_failure = @(sigma1, sigma2, tau12) (sigma1/X_T + sigma2/Y_T + tau12/tau_12f <= 1);
maxstress_failure = @(sigma1, sigma2, tau12) (sigma1/X_T + sigma2/Y_T <= 1);


% Define stress ranges for plotting
sigma_range = linspace(-2*sigma_max, sigma_max, 201);
tau_range = linspace(-2*sigma_max, sigma_max, 201);

% Calculate biaxial stress failure envelopes for Puck and Max Stress criteria
for i = 1:length(sigma_range)
    for j = 1:length(tau_range)
        sigma1 = sigma_range(i);
        tau12 = tau_range(j);
        
        % Calculate sigma2 for [0/90/±45]2s laminate
        a = cosd(45)^2/E1 + sind(45)^2/G12;
        b = sind(45)^2/E1 + cosd(45)^2/G12;
        c = -2*sind(45)*cosd(45)*(1/E1 - 1/G12);
        sigma2 = a*sigma1 + b*tau12 + c*sigma1*tau12;
        
        % Check if stress state satisfies Puck and Max Stress failure criteria
        if puck_failure(sigma1, sigma2, tau12)
            puck_envelope(i,j) = 1;
        else
            puck_envelope(i,j) = 0;
        end
        if maxstress_failure(sigma1, sigma2, tau12)
            maxstress_envelope(i,j) = 1;
        else
            maxstress_envelope(i,j) = 0;
        end
        
        if sigma1 > X_C 
            maxstress_envelope1(i,j) = 1;
        elseif sigma1 < X_T
            maxstress_envelope1(i,j) = 1;
        else
            maxstress_envelope1(i,j) = 0;
        end
        
        if sigma1 > Y_C 
            maxstress_envelope1(i,j) = 1;
        elseif sigma1 < Y_T
            maxstress_envelope1(i,j) = 1;
        else
            maxstress_envelope1(i,j) = 0;
        end



    end
end



% Find FPF and LPF for Puck and Max Stress criteria
for i = 1:length(sigma_range)
    for j = 1:length(tau_range)
        if puck_envelope(i,j) == 1
            puck_LPF_sigma1 = sigma_range(i);
            puck_LPF_tau12 = tau_range(j);
            break
        end
    end
end

% Plot biaxial stress failure envelopes for Puck and Max Stress criteria
figure(1);
hold on;
contour(sigma_range, tau_range, puck_envelope, [0 1], 'b');
contour(sigma_range, tau_range, maxstress_envelope, [0 1], 'r');
xlabel('\sigma_1 (MPa)');
ylabel('\tau_{12} (MPa)');
title('Biaxial Stress Failure Envelopes');
legend('Puck', 'Max Stress');

% Plot FPF and LPF for Puck and Max Stress criteria
scatter(puck_LPF_sigma1, puck_LPF_tau12, 'b', 'filled');
text(puck_LPF_sigma1, puck_LPF_tau12, 'Puck LPF');
legend('Puck', 'Max Stress', 'Puck LPF');

maxstress_LPF_sigma1 = min(sigma_range(maxstress_envelope(end,:)==1));
maxstress_LPF_tau12 = tau_range(end);
scatter(maxstress_LPF_sigma1, maxstress_LPF_tau12, 'r', 'filled');
text(maxstress_LPF_sigma1, maxstress_LPF_tau12, 'Max Stress LPF');
legend('Puck', 'Max Stress', 'Puck LPF', 'Max Stress LPF');



figure(2)
contour(sigma_range, sigma_range, maxstress_envelope1, [0 1], 'r');
xlabel('\sigma_1 (MPa)');
ylabel('\tau_{12} (MPa)');
title('Biaxial Stress Failure Envelopes');
legend('Puck', 'Max Stress');

function maximum_stress_crit(sigma1, sigma2)
end
