function [z_WLC, g_WLC_kbT] = WLC_z_G(f,L)                                   

% extension and internal free energy for stretching the WLC
% input parameters: force ramp (pN) and contour length (nm)
% Artur Kaczmarczyk, June 2018


kbT = 4.114;                                                                 % Boltzmann constant in room temperature (pNnm);
P = 50;                                                                      % persistence length (nm);
S = 1000;                                                                    % stretch modulus (pN)

z_WLC = (L*(1-0.5*sqrt(kbT./(f.*P)))+f./S);                                  % length in nanometers
g_WLC =  z_WLC.*f - L.*f.* (1-sqrt(kbT./(f.*P))+f./(2.*S));                  % internal free energy (pN nm)
g_WLC_kbT = g_WLC./kbT;                                                      % internal free energy (kbT)

end
