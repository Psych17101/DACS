% Rotation Matrix
M_2 = Rotation(theta, M_1);

a = [cos(theta), sin(theta), 0;
    -sin(theta),cos(theta), 0;
    0,0,1];

M_2 = a*M_1;

