%% Main
% Assignement 2023 - DASC - AE4ASM511
% Christophe Hatterer + 
clear; format;
format short g;


%% Initialisation
% Laminate definition (plies of equal thickness)
Nplies = 6;
thetadt = [0 60 -30 30 -60 0]; % ply angles in degrees, from top
thetadb = fliplr(thetadt); % ply angles in degrees, from bottom
h_ply  = 0.003;           % SI units, meters
h      = Nplies * h_ply ;

for i = 1:Nplies
  zbar(i) = - (h + h_ply)/2 + i*h_ply;
end

% Ply engineering properties examples
E1   = 148.e9 ; % Pa
nu12 = .3 ;
E2   = 10.5e9  ; % Pa
G12  = 5.61e9  ; % Pa
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
  [T, invT] = TinvT(theta);
  [Sbar, Qbar] = Tranformed_SCbar(S,Q,T,invT);
  
  A = A + Qbar * h_ply;
  B = B + Qbar * h_ply * zbar(i); 
  D = D + Qbar * (h_ply * zbar(i)^2  + h_ply^3 / 12);
  
  NT = NT + Qbar * h_ply ;
  MT = MT + Qbar * h_ply * zbar(i)  ;
  
  ABD = [A B; A D];
end


%%


