% Equivalent Properties with A matrix
function [E_xA] = EquivPropA(thickness, A_inv)

E_xA = inv(thickness*A_inv(1,1))


