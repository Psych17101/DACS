function [E_x, v_xy, E_y, G_xy, n_xs, n_ys] = elastic_prop_offaxis(E_1,G_12,E_2,v_12,theta)
%
%
%
%
%
%
m = cos(theta);
n = sin(theta);

E_x = 1/(m^4/E_1 + (1/G_12 - 2*v_12/E_1)*m^2*n^2 + n^4/E_2);
E_y = 1/(n^4/E_1 + (1/G_12 - 2*v_12/E_1)*m^2*n^2 + m^4/E_2);

v_xy = E_x*(v_12/E_1*(m^4 + n^4) - (1/E_1 + 1/E_2 - 1/G_12)*m^2*n^2);

G_xy = 1/(2*(2/E_1 + 2/E_2 + 4*v_12/E_1 - 1/G_12)*m^2*n^2 + 1/G_12*(m^4 + n^4));

n_xs = E_x((2/E_1 + 2*v_12/E_1 - 1/G_12)*m^3*n - (2/E_2 + 2*v_12/E_1 - 1/G_12)*m*n^3);
n_ys = E_y((2/E_1 + 2*v_12/E_1 - 1/G_12)*m*n^3 - (2/E_2 + 2*v_12/E_1 - 1/G_12)*m^2*n);
