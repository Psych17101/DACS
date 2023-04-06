function [Sbar, Qbar] = Tranformed_SCbar(S,Q,T,invT)
%
%
%
%
Sbar = invT*S*T;
Qbar = invT*Q*T;
