function [r, v] = integrator_verlet(r, fij_tot)
% intergrator for Newton equations of motion - Verlet Leap-Frog algorithm
global dt m r_old box

%% Verlet algorithm -- update position
r_new = 2*r - r_old + (fij_tot/m)*(dt^2);

% v(t+1)_{x,y} -- velocity calculate
rdiff = (r_new - r_old);
rdiff = rdiff- box*round(rdiff/box);
v = rdiff/(2*dt);

% updating variables
r_old = r;
r = r_new;

end