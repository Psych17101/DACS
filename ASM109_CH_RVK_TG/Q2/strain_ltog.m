function [strain_glo] = strain_ltog(strain_loc,angle)
[T, invT] = TinvT(angle);
% Calculating the global stress vector using the Tinv matrix from local
% stress
strain_glo=invT*strain_loc;