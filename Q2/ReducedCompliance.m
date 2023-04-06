function S_red = ReducedCompliance(E_1,E_2,v_12,G_12)
%ReducedCompliance This function returns the reduced compliance
% matrix for fiber-reinforced materials.
% There are four arguments representing four
% material constants. The size of the reduced
% compliance matrix is 3 x 3.
S_red = [1/E_1 -v_12/E_1 0 ; -v_12/E_1 1/E_2 0;0 0 1/G_12];