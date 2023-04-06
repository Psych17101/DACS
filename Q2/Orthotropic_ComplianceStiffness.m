function [S C] = Orthotropic_ComplianceStiffness(E_1,E_2,E_3,v_12,v_23,v_13, G_12,G_23,G_13)


S_11 = 1/E_1;
S_22 = 1/E_2;
S_33 = 1/E_3;


S_21 = -v_12/E_1;
S_31 = -v_13/E_1;

S_12 = -v_21/E_2;
S_32 = -v_23/E_2;

S_13 = -v_31/E_3;
S_23 = -v_32/E_3;

S_44 = 1/G_23;
S_55 = 1/G_13;

S_66 = 1/G_12;


S = [S_11 S_12 S_13 0 0 0;
    S_21 S_22 S_23 0 0 0; 
    S_31 S_32 S_33 0 0 0;
    0 0 0 S_44 0 0;
    0 0 0 0 S_55 0;
    0 0 0 0 0 S_66];
C = inv(S);


