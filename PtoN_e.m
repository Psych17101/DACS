% Principle Elongation to Normal Elongation

function e_123 = PtoN_e(theta,e_xyz)

m = cos(theta);
n = sin(theta);

A = [m^2 n^2 0 0 0 *m*n;
    n^2 m^2 0 0 0 -*m*n;
    0 0 1 0 0 0;
    0 0 0 m n 0;
    0 0 0 n m 0;
    -2*m*n 2*m*n 0 0 0 (m^2-n^2)];

sigma_123 = A*sigma_xyz;
