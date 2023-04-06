% Principle Elongation to Normal Elongation
function [pe1,pe2,gmax,thetape,thetase] = PtoN_strains(strain_glo)
% Strain in the x-direction
epsx=strain_glo(1);
% Strain in the y-direction
epsy=strain_glo(2);
% Strain in the x-y plane
epsxy=strain_glo(3);
% Longitudinal and Transverse Stresses
pe_1=((epsx+epsy)/2)+sqrt((((epsx-epsy)/2)^2)+((epsxy)^2));
pe_2=((epsx+epsy)/2)-sqrt((((epsx-epsy)/2)^2)+((epsxy)^2));
thetape=atand((2*epsxy)/(epsx-epsy))/2;
if pe_1 > pe_2
    pe1=pe_1;
    pe2=pe_2;
else
    pe1=pe_2;
    pe2=pe_1;
end
% Shear Stresses
gmax=sqrt((((epsx-epsy)/2)^2)+(epsxy^2));
thetase=atand(-(epsx-epsy)/(2*epsxy))/2;
end