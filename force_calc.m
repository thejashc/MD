function [fij] = force_calc(rij, fij)

global epsilon sigma rc2 N pot_choice

% final shortest distances between the particles
rij_fin = sqrt(rij(:,:,1).^2 + rij(:,:,2).^2);

if (pot_choice==0)
    %% force calculation -- Quadratic potential
    dphi_drij = -2*epsilon*(sigma - rij_fin);
    fij(:,:,1) = -dphi_drij.*(rij(:,:,1)./rij_fin);
    fij(:,:,2) = -dphi_drij.*(rij(:,:,2)./rij_fin);
    
elseif (pot_choice==1)
    
    %% force calculation -- Lennard-Jones fluid
    dphi_drij = (24*epsilon*(rij_fin).^(-1)).*( -2*(rij_fin/sigma).^(-12) + (rij_fin/sigma).^(-6));
    fij(:,:,1) = -dphi_drij.*(rij(:,:,1)./rij_fin);
    fij(:,:,2) = -dphi_drij.*(rij(:,:,2)./rij_fin);
end

% force threshold in distance
fij(isnan(fij))=0;                                       % set values above threshold to 0
ind = find(rij_fin >= rc2);                               % indices for rij_fin > rc2
ind2 = ind + N^2;                                        % generate 2nd set of indices
fij(ind) = 0;                                            % threshold force for x,y components
fij(ind2)=0;

end