function [rij] = rij_calc(r, rij)

global N box

% initialize positions
[i,j] = meshgrid(1:N, 1:N);
rij(:,:,1) = reshape(r(i,1) - r(j,1), N,N);
rij(:,:,2) = reshape(r(i,2) - r(j,2), N,N);

% Periodic boundary condition -- shortest x and y distance between particles
rij = rij- box*round(rij/box);

end

