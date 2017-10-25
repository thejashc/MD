function [r, w, v] = integrator_NoseHoover(r, w, fij_tot)
% intergrator for Newton equations of motion - Verlet Leap-Frog algorithm
% implementing the Nose-Hoover dynamics 
% first updating v(t-dt/2) -> v(t+dt/2) 
% use v(t+dt/2) to update r(t-dt) -> r(t), eta(t-dt) -> eta(t)

global dt eta Q dof T_set

% Verlet Leap-frog scheme
% storing v(t - dt/2) in wold
wold = w;

% obtaining w or v(t + dt/2) from r(t) (fij_tot), eta(t) and v(t - dt/2) 
% using implicit formula: report has the derivation
w = (wold*(1 - 0.5*eta*dt) + dt*fij_tot)/(1 + 0.5*eta*dt);

% updating position using v(t+dt/2)
r = r + w*dt;

% updating viscosity using v(t+dt/2)
w2sq = sum(sum(w.^2));
eta = eta + (dt/Q)*(w2sq - dof*T_set);

% calculating v(t) from v(t - dt/2) and v(t + dt/2) 
v = (w + wold)/2.0;

end
