%% Main
% Assignement 2023 - DASC - AE4ASM511
% Christophe Hatterer + Rasmus van Kerkvoorde
clear; format;
format short g;


%% Initialisation
% Laminate definition (plies of equal thickness)

thetadt = [0 +45 -45 90 90 -45 +45 0]; % ply angles in degrees, from top??? Not bottom?
Nplies = length(thetadt);
thetadb = fliplr(thetadt); % ply angles in degrees, from bottom
h_ply  = 0.125*10^(-3);           % SI units, meters
h      = Nplies * h_ply ;

z = 0:h_ply:h;

for i = 1:Nplies
  zbar(i) = - (h + h_ply)/2 + i*h_ply;
end

% Ply engineering properties examples
E1   = 140.e9 ; % Pa
nu12 = .3 ;
E2   = 10.e9  ; % Pa
G12  = 5.e9  ; % Pa
nu21 = nu12 * E2 / E1 ;
X_T = 150.e6;
X_C = 120.e6;
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


for i = 1:Nplies
  % For each ply we calculate the ABD Matrix
  theta  = thetadb(i)*pi/180;  % ply i angle in radians, from bottom
  [T, invT] = TinvT(theta); %% there may be a problem with the TinvT and Transformed_SCbar
  [Sbar, Qbar] = Tranformed_SCbar(S,Q,T,invT);
  
  A = A + Qbar * (z(i+1)-z(i)) ; %N/m, right dimensions?
  B = B + (1/2)*Qbar * (z(i+1)^2-z(i)^2); %N 
  D = D + (1/3)*Qbar * (z(i+1)^3-z(i)^3); %Nm
  
  A_inv = inv(A);
  B_inv = inv(B);
  D_inv = inv(D);

 

  NT = NT + Qbar * (z(i+1)-z(i));
  MT = MT + (1/2)*Qbar * (z(i+1)^2-z(i)^2) ;
  
  %changes on every iteration, and isnt it supposed to be [A B;B D]
  %there is only one ABD matrix, meaning line 63 should be outside of the
  %loop
  ABD = [A B; A D];
end

% Exercise done is class
E_x = EquivPropA(h,A_inv); 
N_initial  = [0 20000]; %N/cm
iter1   = 0;
k3      = 0;



while N_initial(2)-N_initial(1)>0.01
    N_x = mean(N_initial);
    N_t = [N_x 0 0];

    eps_glo = A_inv*N_t';
    FI = 0;

    for i = 1:Nplies
        % For each ply we calculate the ABD Matrix
        theta  = thetadb(i)*pi/180;  % ply i angle in radians, from bottom
        [eps_loc] = strain_gtol(eps_glo,theta);
        disp(eps_loc)
        sigma_loc = Q*eps_loc';
        disp(sigma_loc)

        %Failure test
        [F1,F2,F3] = Failure1(sigma_loc,X_T,X_C,Y_T,Y_C,S_t);
        Maxi_F = max([F1,F2,F3]);
        disp(Maxi_F)
        if Maxi_F > FI
            FI = Maxi_F;
        end    
    end
        N_xnew = N_t(1)/FI;
        if N_xnew > N_initial(1) && N_xnew <N_initial(2)
            N_initial(1) = N_xnew;
        elseif N_xnew > N_initial(2)
            N_initial(2) = N_xnew;
        end

end
    
% end



%% Functions

function [T, invT]= TinvT(theta)
%T This function returns the transformation matrix T
% given the orientation angle "theta".
% There is only one argument representing "theta"
% The size of the matrix is 3 x 3.
% The angle "theta" must be given in degrees.
m = cos(theta);
n = sin(theta);
T = [m*m n*n 2*m*n ; n*n m*m -2*m*n ; -m*n m*n m*m-n*n];
invT = [m*m n*n -2*m*n ; n*n m*m 2*m*n ; m*n -m*n m*m-n*n];
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
[T, invT] = TinvT(angle);
% Calculating the global stress vector using the Tinv matrix from local
% stress
strain_glo=invT*strain_loc;
end

function [strain_loc] = strain_gtol(strain_glo,angle)
[T, invT] = TinvT(angle);
% Calculating the global stress vector using the Tinv matrix from local
% stress
strain_loc=strain_glo'*T;
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



