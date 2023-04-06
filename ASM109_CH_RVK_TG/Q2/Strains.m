function y = Strains(eps_xo,eps_yo,gam_xyo,kap_xo,kap_yo,kap_xyo,z)
% Strains This function returns the strain vector at any point P
% along the normal line at distance z from point Po which
% lies on the reference surface. There are seven input
% arguments for this function - namely the three strains
% and three curvatures at point Po and the distance z.
% The size of the strain vector is 3 x 1.
epsilonx = eps_xo+z* kap_xo;
epsilony = eps_yo+z* kap_yo;
gammaxy = gam_xyo+z* kap_xyo;
y = [epsilonx ; epsilony ; gammaxy];