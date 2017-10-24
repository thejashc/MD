clearvars;
clc;

global N epsilon sigma box steps dt m rc2 kB d pot_choice r_old phi_lj_rc2

box = 10;
N = box*box;
steps = 1e5;
dt = 1e-4;
m = 1;
kB = 1;
d = 2;
dr = 0.25;
trial = 1;
MSD_lim = 5e4;
rad = 0.5:0.25:10;
pot_choice = 1;

% choice of potential
if (pot_choice==0)
    
    %% quadratic-soft potential
    epsilon = 100;
    sigma = 3.0;
    rc2 = 3.0;
elseif (pot_choice ==1)
    
    %% Lennard Jones
    epsilon = 1.0;
    sigma = 1.0;
    rc2 = 2.5;
    
    phi_lj_rc2 = 4*epsilon*( (sigma/rc2)^12 - (sigma/rc2)^6 );
end


for t=1:trial
    
    % declaring variables
    r = box*(rand(N,2) - 0.5);
    [posx,posy]=meshgrid(1:box, 1:box);
    rinit = [reshape(posx, N,1), reshape(posy, N,1)];
    xcom = (1/N)*sum(rinit(:,1));
    ycom = (1/N)*sum(rinit(:,2));
    r = rinit - [xcom,ycom];
    
    v = zeros(N,2);
    fij = zeros(N, N, 2);
    rij = zeros(N, N, 2);
    T = zeros(steps,1);
    P = zeros(steps,1);
    N_count = zeros(length(rad),1);
    phi_tot = zeros(steps, 1);
    kin_tot = zeros(steps, 1);
    e_tot = zeros(steps,1);
    MSD = zeros(MSD_lim,1);
    n = zeros(N,2);
    count = 0;
    msd_count =0;
    k=0;
    vdist = [];
    rdist = [];
    
    % calculating rij matrix
    rij = rij_calc(r, rij);
    
    for step=1:steps
        %% calculation fij between particles
        fij = force_calc(rij, fij);
        
        % total x and y component of force calculation
        fij_tot = [ sum(fij(:,:,1),1)', sum(fij(:,:,2),1)' ];
        
        if (step == 1)
            w = rand(N,2);
            w = w - mean(w);
            
            % v = importdata('vel_start_up.dat');
            v = rand(N,2);
            v = v - mean(v);
            r_old = r - v*dt;
        end
        
        %% equation of motion integrator
        [r, w, v] = integrator_verlet_lf(r, w, fij_tot);
        % [r, v] = integrator_verlet(r, fij_tot);
        % [r, v] = integrator_euler(r, v, fij_tot);
        
        % updating the n vector and calculating rij matrix
        indpx = find( r(:,1) > box/2 );
        indnx = find( r(:,1) < -box/2 );
        indpy = find( r(:,2) > box/2 );
        indny = find( r(:,2) < -box/2 );
        n(indpx,1) = n(indpx,1) + 1;
        n(indnx,1) = n(indnx,1) - 1;
        n(indpy,2) = n(indpy,2) + 1;
        n(indny,2) = n(indny,2) - 1;
        
        %% checks
        %if ( ~isempty(indpx) || ~isempty(indnx) || ~isempty(indpy) || ~isempty(indny))
        %    flag = 1;
        %end
        mom = mean(v);
        if (mom(1) > 1 || mom(2) > 1 || step==17170)
            flag = 1;
        end
        
        %% periodic bc's
        r = mod(r+ box/2, box) - box/2;
        rij = rij_calc(r, rij);
        
        % measuring conserved quantities
        [phi_tot(step), kin_tot(step), e_tot(step)] = conserve(rij, v);
        
        %% measuring pressure and temperature
        % temperature
        dof = d*(N-1) - 1;
        T(step) = 2*kin_tot(step)/(dof*kB);
        
        % % pressure
        rF = rij(:,:,1).*fij(:,:,1) + rij(:,:,2).*fij(:,:,2);
        P(step) = (1/box^2)*( N*kB*T(step) +  (0.5/d)*sum(sum(rF,2),1));
        
        %% velocity and g(r) calculation
        if ((step > 2e4) && (mod(step, 100) == 0) )
            % velocity distribution
            vdist = [vdist; v];
            rdist = [rdist; r];
            
            % g(r) calculation
            % final shortest distances between the particles
            rij_fin = sqrt(rij(:,:,1).^2 + rij(:,:,2).^2);
            
            for i=1:N
                for j=1:length(rad)
                    N_count(j) = N_count(j) + length( find( rij_fin(:,i) > rad(j) & rij_fin(:,i) < rad(j)+dr) );
                end
                
                % count variable
                count = count + 1;
            end
            
        end
        
        %% MSD calculation
        if( mod(step,MSD_lim) == 0 )
            rinit = r+n*box;
            k= 1;
            
            % msd_count
            msd_count = msd_count + 1;
        end
        
        if ( (step > MSD_lim) && mod(step - k, MSD_lim) == 0 && (k <= MSD_lim) )
            MSD(k) =  MSD(k) + sum(sum(((r+ n*box) - rinit ).^2));
            
            % update k
            k = k+1;
        end
        
        %visualization
        % visualize(r);
        % pause(0.01)
    end
    % verifying against ideal-gas law
    P_mean = mean(P(end-3000:end));
    T_mean = mean(T(end-3000:end));
    
    id_ratio(t)= P_mean/((N/box^2)*(kB*T_mean));
    
    %% Energy conservation
    figure(1)
    semilogx(1:steps, phi_tot, '.-')
    hold on
    semilogx(1:steps, kin_tot,'.-')
    semilogx(1:steps, e_tot,'.-')
    
    %     %% Velocity distribution
    [y,x] = hist(vdist(:,1),50);
    y = y/trapz(x,y);
    figure(2);
    plot(x, y, 'b-o')
    s = -60:1e-3:60;
    gauss = (1/sqrt(2*pi*kB*T_mean))*exp(-s.^2 / (2*kB*T_mean));
    hold on;
    plot(s, gauss, 'r-')
    
    
    %% plotting g(r)
    figure(3);
    rho_homo = N/(box^2);
    g = N_count./(rho_homo*2*pi*rad'*dr*count);
    
    plot(rad, g, 'k-')
    
    %% plotting the MSD
    figure(4);
    MSD_t = 1:MSD_lim;
    MSD_y = MSD/(N*(msd_count-1));
    
    xfit= log(MSD_t(end-50:end-1));
    yfit = log(sqrt(MSD_y(end-50:end-1)));
    
    plot(xfit, yfit, 'r^-')
    
    t
end

% writing data
% csvwrite(['./data/dt_',num2str(dt),'_phi_tot.dat'],[(1:steps)' phi_tot]);
% csvwrite(['./data/dt_',num2str(dt),'_kin_tot.dat'],[(1:steps)' kin_tot]);
% csvwrite(['./data/dt_',num2str(dt),'_e_tot.dat'],[(1:steps)' e_tot]);
% csvwrite('./data/N_200_box_50_sig_1_0_vdist_theory_compare.dat',[x' y']);
% csvwrite('./data/N_200_box_50_sig_1_0_energy_balance.dat',[(1:1e5)' phi_tot kin_tot e_tot]);
