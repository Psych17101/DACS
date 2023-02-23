% Rotation Matrix
function T = Rotation(theta)
%T This function returns the transformation matrix T
% given the orientation angle "theta".
% There is only one argument representing "theta"
% The size of the matrix is 3 x 3.
% The angle "theta" must be given in degrees.
m = cos(theta);
n = sin(theta);
T = [m*m n*n 2*m*n ; n*n m*m -2*m*n ; -m*n m*n m*m-n*n];
