function [r, v] = integrator_euler(r, v, fij_tot)
% intergrator for Newton equations of motion - Verlet Leap-Frog algorithm
global dt m 

% update velocities
v = v + (fij_tot/m)*dt;

% update positions using first order 
r = r + v*dt;

end