function [stress_glo] = stress_ltog(stress_loc,angle)
[T, invT] = TinvT(angle);
% Calculating the global stress vector using the Tinv matrix from local
% stress
stress_glo=invT*stress_loc;