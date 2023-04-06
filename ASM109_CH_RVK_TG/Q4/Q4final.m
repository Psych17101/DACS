clear; format;
close all;
format short g;

%% Initialisationation
layup = [0 90 +45 -45 -45 +45 90 0 0 90 +45 -45 -45 +45 90 0];
Nplies = length(layup);
thetadb = fliplr(layup); % ply angles in degrees, from bottom;
h_ply  = 0.12784*10^(-3);           % SI units, m
h      = Nplies * h_ply ;

z = -h/2:h_ply:h/2;

for i = 1:Nplies
  zbar(i) = - (h + h_ply)/2 + i*h_ply; % mm
end

% material properties (mean and standard deviation) in Pa
mean_E1 = 165.22e9; std_E1 = 16.922e9;
mean_E2 = 8.444e9; std_E2 = 1.077e9;
mean_G12 = 6.412e9; std_G12 = 1.109e9;
mean_v12 = 0.3532; std_v12 = 0.1809;
mean_X_T = 1924e6; std_X_T=108.7*1e6;% Fiber
X_C = 1480e6; %
mean_Y_T = 107.2e6; std_Y_T=9.351e6;% Matrix
Y_C = 220e6;
mean_S12 = 152.4e6; std_S12=1.784e6;


% load and orientation
N = 400*1e3; % resultant force in N/m

theta = 30; % angle of resultant force with respect to X-axis in degrees

% Monte Carlo simulation parameters

F = [N*cosd(theta); N*sind(theta); 0]; % load vector N/m
error=1;
pflist=[];
hh = animatedline('Color','r');
axis([0,10000,0,0.3])
FN=1;
FNlast=1;
n_simulations=0;
fel=[];
while error>1e-8 && n_simulations<10000 %run until converged
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
        S12 = normrnd(mean_S12, std_S12, 1);
        % calculate the transformed stiffness matrix for each sample

        %strain_samples = zeros(3);
        %stress_samples = zeros(3);
   
        A = zeros(3,3);
        B = zeros(3,3);
        D = zeros(3,3);
        [S, Q] = ReducedComplianceStiffness(E1_samples,E2_samples,v12_samples,G12_samples);
        [moduli] = [E1_samples,E2_samples,v12_samples,G12_samples]; 
        for i = 1:Nplies
            % For each ply we calculate the ABD Matrix
            [Qbar,Sbar] = QbarandSbar(thetadb(i),moduli);

            A = A + Qbar * (z(i+1)-z(i)) ; %N/mm
            B = B + (1/2)*Qbar * (z(i+1)^2-z(i)^2); % N 
            D = D + (1/3)*Qbar * (z(i+1)^3-z(i)^3); %N mm
            ABD = [A B; B D];
        end
        

       
        strain_samples = inv(A)*F;%invA_samples*F; ;
        stress_samples = Qbar*strain_samples;
        for l = 1:Nplies
            % Calculations of Strains and Stresses
             % global sigmaxx % Pa
            [eps_loc] = strain_gtol(strain_samples(:),thetadb(l)); % ply i angle in radians, from bottom
            sigma_loc = stress_gtol(stress_samples, thetadb(l));


            % Calculations of Puck criterion for each simulation run,
            % each random variable chosen and each ply 
            %fe = puckfailure(sigma_loc,E1_samples,E_x,v12_samples,X_T,X_C);  % (num_samples,layer)
            [fe1, fe2] = PuckCriterion(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,S12,v12_samples,E1_samples,'c','n');
            % Calculation of Ply failure 
            % Calculation of max failure index under puck criterion
                % for each simulation and each random sample
            %if fe > max_fe
            %    max_fe = fe; % failure index under puck criterion
            %    %numberply(j,k,l) = l;
            %end
            nfe=max(fe1, fe2);
            if nfe>fe
                fe=nfe;
            end
        
        end 
        
        if fe<1
            Truth=0;
            
            %disp("JIPIIE")
        end
    fel=[fel, fe];
    if length(fel)>1000
        Truth=0;
        mean(fel)
        s=fel;
        disp("done")
        
        fel=[];
    end
    fe=0;
    n=n+1;
    end
    pflist=[pflist, 1/n];%P_f^R
    % is the failure index bigger than 1?
    %FNlast2=FNlast;
    %FNlast=FN;
    FN=sum(pflist)/length(pflist);
    
    %plot(length(FN), FN); hold on;
    %ROC=abs(FNlast2-FNlast)/abs(FNlast-FN);
    addpoints(hh,length(pflist),FN);
    
    addpoints(hh, length(pflist), ROC); 
    drawnow;
    %error=abs(FN-FNlast);
    n_simulations=n_simulations+1;
    
    
end

% calculate the probability of failure for each load and orientation
%P_f_convergence(k) = sum(ply_failure) / num_simulations;



%%
S100=[0.24388, 0.2423, 0.24146, 0.24447, 0.24782, 0.24251, 0.24315, 0.23599, 0.24232, 0.24446
];
meanS100=mean(S100); %0.24284
stdS100=std(S100); %0.0029955
x=meanS100-4*stdS100:0.0001:meanS100+4*stdS100;
plot(x, normpdf(x, meanS100, stdS100));xlabel("Probability of failure"); grid on;

t_stat = tinv(0.995, length(S100)-1);

% Calculate the margin of error
margin_error = t_stat * stdS100 / sqrt(length(S100));

% Calculate the lower and upper bounds of the confidence interval
lower_bound = meanS100 - margin_error; %0.23976
upper_bound = meanS100 + margin_error; %0.24591
%% 
S400=[3.5614, 3.5548, 3.5670, 3.5630, 3.5573, 3.5551, 3.559, 3.5565, 3.5618, 3.5565];
meanS400=mean(s); %0.24284
stdS400=std(s); %0.0029955
x=meanS400-4*stdS400:0.0001:meanS400+4*stdS400;
plot(x, normpdf(x, meanS400, stdS400));xlabel("Failure index"); grid on;
prob = normcdf(1, meanS400, stdS400);
%%
probl=[7.4404e-16, 1.1047e-15, 8.4353e-16, 1.0793e-15, 6.6061e-16, 7.8039e-16, 1.0115e-15,1.1453e-15, 1.0931e-15,1.2127e-15];
mean1=mean(probl);
std1=std(probl);


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


                    %PuckCriterion(sigma_loc(1),sigma_loc(2),sigma_loc(3),X_T,X_C,Y_T,Y_C,S12,v12_samples,E1_samples,'c','n')
function [FF,IFF] = PuckCriterion(sigma1,sigma2,sigma3,X_T,X_C,Y_T,Y_C,G12_t,nu12,E1,fiber,print)
FF =0;
IFF=0;
sigma3 = abs(sigma3);

% Determination of which fiber used
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
%IFF
end