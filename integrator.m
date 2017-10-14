function [r, w, v] = integrator(r, w, fij_tot)
% intergrator for Newton equations of motion - Verlet Leap-Frog algorithm
global dt m

% saving variables
% temp=r;

% % Verlet algorithm -- update position
% r(:,1) = 2*r(:,1) - rpre(:,1) + (fij_tot(:,1)/m)*(dt^2);
% r(:,2) = 2*r(:,2) - rpre(:,2) + (fij_tot(:,2)/m)*(dt^2);
% 
% % v(t+1)_{x,y} -- velocity calculate
% v(:,1) = (r(:,1) - rpre(:,1))/(2*dt);
% v(:,2) = (r(:,2) - rpre(:,2))/(2*dt);

% v(:,1) = v(:,1) + fij_tot(:,1)*(dt/m);      % velocity x-component
% v(:,2) = v(:,2) + fij_tot(:,2)*(dt/m);      % velocity y-component
% 
% % r(t+1)_{x,y} -- position update
% r(:,1) = r(:,1) + v(:,1)*dt;                % x-component of posn
% r(:,2) = r(:,2) + v(:,2)*dt;                % x-component of posn

% updating rpre
% rpre = temp;

%% Verlet Leap-frog scheme
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

