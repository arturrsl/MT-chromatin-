function [z_final, f] = fit_FE_fiber(f,CL,Nnuc,NRL,kf,degeneracy,dG1,dG2)  

% simulating chromatin fiber force-extension curves

% based on: Meng, Andresen, van Noort - NAR, 2015
% modification with respect to the paper: transition from unwrapped to
% extended state is described only with WLC function

% input parameters: 
% force ramp (pN),contour length (nm), number of
% nucleosomes, Nucleosome Repeat Length (bp), fiber stiffness (pN/nm), degeneracy, 
% free energy of unstacking (dG1), free energy of the intermediate transition (dG2)

% Artur Kaczmarczyk, June 2018

%%

%f = (0.03:0.01:7);
kbT = 4.114;                                     % Boltzmann constant in room temperature (pNnm);
%Nnuc = 15;                                      % number of assembled nucleosomes
%NRL = 167;                                      % Nucleosome Repeat Length (bp)
%CL = 4535;                                      % length of the DNA template used for reconstitution(bp)
Lextended = 5;                                   % additional gain in extension after unstacking and 1st nucleosomal turn unwrapping (nm/nuc)
%kf = 1;                                         % stiffness of a folded chromatin fiber per nucleosome (pN/nm)
Lwrap = 89;                                      % length of DNA unwrapped from the first nucleosomal turn (bp)
%degeneracy = 0;                                

%dG1 = 22;                                       % free energy of unstacking (dG1)
%dG2 = 11;                                       % free energy of the intermediate transition (dG2)
dG3 = 100;                                       % high force transition (not in equlibrium, value just for plotting purposes)


L_tether = CL - Nnuc * NRL;                      % length of bare DNA handles (bp)

%% help plots with transition borders

[zet_singlywrapped] = WLC_z_G(f,(CL - Nnuc.*Lwrap).*0.34 );             
[zet_extended] = WLC_z_G(f,(CL - Nnuc.*Lwrap).*0.34  + Nnuc.*Lextended);
[zet_unwrapped] = WLC_z_G(f,CL .* 0.34);

plot(zet_singlywrapped./1000,f,':');
hold on;
plot(zet_extended./1000,f,':');
plot(zet_unwrapped./1000,f,':');

%% tether length at different states (per nucleosome)

[z_wrapped, G_wrapped] = WLC_z_G(f,L_tether .* 0.34);                       % total DNA handles (nm)

[z_fiber, G_fiber] = fiber_z_G(f,kf);                                       % fiber extension and internal free energy per nucleosome [nm]

[z_singlewrap, g_singlewrap] = WLC_z_G(f,(NRL - Lwrap) .* 0.34);            % DNA extension and internal free energy of unwrapped base pairs; per nucleosome [nm]
G_singlewrap = g_singlewrap + dG1 ;                                         % correction with free energy G1 that overcomes the energy barrier

[z_extended, g_extended] = WLC_z_G(f,(NRL - Lwrap) .* 0.34 + Lextended);    % DNA extension and internal free energy of unwrapped base pairs; per nucleosome [nm]
%G_extended = G_singlewrap +  dG2 - f.*Lextended./kbT;                      % correction with free energy G2 that overcomes the energy barrier and the work done due to the transition (in other transitions, the work done is included in the WLC/fiber functions)
G_extended = g_extended +  +dG1 + dG2;                                      % correction with free energy G2 that overcomes the energy barrier

[z_unwrapped, g_unwrapped] = WLC_z_G(f,NRL .* 0.34);                        % DNA extension and internal free energy of unwrapped base pairs; per nucleosome [nm]
G_unwrapped = g_unwrapped + dG1 + dG2 + dG3;                                % correction with free energy G3 that overcomes the energy barrier

%% calculation of all possible states that nucleosome in an array can form, its extension and free energy (!!!work is substracted -z * F !!!)
%  rows: states, columns: elements of the force ramp

D = [];
state = [];
s = 1;

for i = 1:(Nnuc+1)
   
    for j = 1: (Nnuc - i +3)  
        
        for k = 1: (Nnuc - i - j +3)
            
               state(s,:) = [i-1,j-1,k-1,Nnuc - (i-1) - (j-1) - (k-1)];
               D(s) = factorial(Nnuc)./(factorial(i-1).*factorial(Nnuc-(i-1))) .* (factorial(Nnuc-(i-1))./(factorial(j-1).*factorial(Nnuc-(i-1)-(j-1)))) .* (factorial(Nnuc-(i-1)-(j-1))./(factorial(k-1) .*factorial(Nnuc-(i-1)-(j-1)-(k-1)))); % this particular state exists in D combinations
               D(s) = 1+(D(s)-1) .* degeneracy;                             % a trick to switch off the degeneracy if not needed                      
            
               D_array(s,:) = repmat(D(s),1,length(f));                     % copying the degeneracy to all columns                 
               
               z_tot(s,:) = z_wrapped + state(s,1) .* z_fiber + state(s,2) .* z_singlewrap + state(s,3) .* z_extended + state(s,4) .* z_unwrapped;          % extension of the states  multiplied by the amount of elements in this particular state
            
               G_tot(s,:) = ((G_wrapped + state(s,1) .* G_fiber + state(s,2) .* G_singlewrap + state(s,3) .* G_extended + state(s,4) .* G_unwrapped)) - f.*z_tot(s,:)./kbT;     % internal energy + energy barriers + work done by the tether (-z * F)
               
    
               s = s+1;
               
               
        end
        
    end
    
end

%% Boltzmann weighing - the state with the lowest total free energy dominates in total extension

for i = 1:length(f)
    
    G_total(:,i) = G_tot(:,i) - min(G_tot(:,i));                            % offseting the absolute values of energies to avoid large negative numbers that go to exponent

end
    

Bolz_top = z_tot .* D_array .* exp (-G_total);                              
Bolz_bot = D_array .* exp (-G_total);

[row col] = size(z_tot);

for i = 1:col
    z_final(i) = (sum(Bolz_top(:,i))./sum(Bolz_bot(:,i))./1000);
   
end

plot(z_final, f)
hold on;
