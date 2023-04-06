function [S, Q] = ReducedComplianceStiffness(E_1,E_2,v_12,G_12)
%ReducedStiffness This function returns the reduced stiffness
% matrix for fiber-reinforced materials.
% There are four arguments representing four
% material constants. The size of the reduced
% stiffness and compliance matrix is 3 x 3.
v_21 = v_12*E_2/E_1;
S = [1/E_1 -v_12/E_1 0 ;
    -v_12/E_1 1/E_2 0;
    0 0 1/G_12];
Q = [E_1/(1-v_12*v_21) v_12*E_2/(1-v_12*v_21) 0 ;
    v_12*E_2/(1-v_12*v_21) E_2/(1-v_12*v_21) 0 ; 
    0 0 G_12];

