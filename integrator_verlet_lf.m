function [r, w, v] = integrator_verlet_lf(r, w, fij_tot)
% intergrator for Newton equations of motion - Verlet Leap-Frog algorithm
global dt m

% Verlet Leap-frog scheme
wold = w;

% updating w
w(:,1) = w(:,1) + fij_tot(:,1)*(dt/m);
w(:,2) = w(:,2) + fij_tot(:,2)*(dt/m);

% updating position
r(:,1) = r(:,1) + w(:,1)*dt;
r(:,2) = r(:,2) + w(:,2)*dt;

% calculating velocity
v = (w + wold)/2.0;

end