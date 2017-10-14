function [phi_tot, kin_tot, e_tot] = conserve(rij, v)
% calculating the conserved quantities

global epsilon sigma m rc2 pot_choice

% final shortest distances between the particles
rij_fin = sqrt(rij(:,:,1).^2 + rij(:,:,2).^2);
ind = (rij_fin >= rc2);                                                      % indices for rij_fin > rc2

% calculating the total potential
if (pot_choice ==0)
    
    %% quadratic potential
    phi_rij = epsilon*(sigma - rij_fin).^2;                                     % calculating the potential for particles (i,j)
    
elseif (pot_choice == 1)
    
    %% lennard jones fluid
    phi_rij = 4*epsilon*( (rij_fin/sigma).^(-12) - (rij_fin/sigma).^(-6));  
    phi_rij(isnan(phi_rij))=0;                                       % set values above threshold to 0
    
end

phi_rij(ind) = 0;
phi_tot = 0.5*(sum(sum(phi_rij,1),2) - trace(phi_rij));                     % calculating the total potential energy

% calculating the total kinetic energy
kin_tot = 0.5*m*sum(sum(v.^2,2),1);

% total energy
e_tot = phi_tot + kin_tot;

end