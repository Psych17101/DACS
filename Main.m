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
N_initial  = [0 100000]; %N/m
iter1   = 0;
k3      = 0;

while I(2)-I(1)>0.01%im guessing this is error
    mean = mean(N_initial);

end


disp(A)



%% Functions




