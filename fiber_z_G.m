function [z_fiber, g_fiber_kbT] = fiber_z_G(f,kf)       

% extension and internal free energy for stretching the Hookean spring 
% input parameters: force ramp (pN) and stretching stiffness (pN/nm)
% Artur Kaczmarczyk, June 2018


kbT = 4.114;                                                                 % Boltzmann constant in room temperature (pNnm);
z0 = 1.7;                                                                    % spring initial length (chromatin fiber length per nucleosome) in the absence of force                                                                                                                                        % z0 - fiber length per nucleosome

z_fiber = f ./ kf + z0;                                                      % extension in nanometers
g_fiber =   f .^ 2 ./ (2 .* kf) - f .* z0;                                   % internal free energy (pNnm)
g_fiber_kbT = g_fiber./kbT;                                                  % internal free energy (kT)

end
