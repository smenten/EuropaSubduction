


%% Friction equation (temperature change) **WORKING**

k = 2.3  ;              % thermal conductivity of water ice 
rho_i = 917 ;           % density of ice in kg/m^3
C = 2050  ;             % specific heat capacity of water ice
u = 1e-15:1e-15:6.342e-9;              % slab velocity (m/s) (between 0 and 200 mm/yr)
phi = 0;                % angle between direction of velocity and slab subduction
theta = 15 ;           % angle of subduction 
miu = 0.7 ;            % coefficient of friction (ice sliding on ice)
d_v = 3000 ;            % depth to melt
g = 1.315 ;             % gravity of Europa in m/s^2 
kappa = k/(rho_i*C) ;       % thermal diffusivity of water ice


tau = miu.*rho_i.*g.*d_v     ;                 % shear stress

a = ((2.*tau)./k) ;
b = (kappa.*d_v.*u) ;
c = (pi.*cosd(phi).*sind(theta)) ;

d_T = a.*((b./c).^(0.5)) ;

u_n = u.*3.154e10 ; % velocity converted from m/s to mm/yr

plot (u_n, d_T, 'LineWidth', 2)
hold on
% new angle
theta = 30;

a = ((2.*tau)./k) ;
b = (kappa.*d_v.*u) ;
c = (pi.*cosd(phi).*sind(theta)) ;

d_T = a.*((b./c).^(0.5)) ;
u_n = u.*3.154e10 ;
plot (u_n, d_T , 'LineWidth', 2)
% new angle


theta = 45;

a = ((2.*tau)./k) ;
b = (kappa.*d_v.*u) ;
c = (pi.*cosd(phi).*sind(theta)) ;

d_T = a.*((b./c).^(0.5)) ;
u_n = u.*3.154e10 ;
plot (u_n, d_T , 'LineWidth', 2)

hold off

xline(40, '--', 'LineWidth', 2)
xlabel('Plate Velocity (mm/yr)')
ylabel('Change in Temperature (K) (at 10 km depth)')
xlim ([0 200]);
ylim ([0 60]);
set(gca,'FontSize',15)
legend ('15{\circ}', '30{\circ}', '45{\circ}')


% %% For-Loop Time Dependent change in T, not working
% 
% % Constants
% Temp =  100;            % starting temperature of subducting slab at surface
% T = Temp ;
% k_t = 651/T;            %temperature dependent thermal conductivity of ice (Johnson et al 2017)
% c_t = 1925*(T/250) ;             % temperature dependent thermal conductivity of ice (Johnson et al 2017)
% g = 1.315 ;             % gravity of Europa in m/s^2
% u = 1.268e-9;              % slab velocity (m/s)  (40 mm/yr)
% C = 2050  ; 
%             % from temperature of conductive portion of the lid in Kelvin, is linear with depth
% u_n = u.*3.154e10 ; % velocity converted from m/s to mm/yr  
% rho_i = 917 ;           % density of ice in kg/m^3
% kappa = k_t/(rho_i*C) ;       % thermal diffusivity of water ice
% 
% 
% %okay: so every year = 40 mm change in depth on the hypotenuse
% %
% % coefficient of friction changing as well 
% 
% 
% sub_A = 30;           % subduction angle
% hypot = 0.04 ;      % 40 mm, aka moved 40 mm along the hypotenuse to start
% 
% 
% 
% 
% 
% for z = 1:200    % time in years
% 
% i=2;
% 
% k_t =  2.3; %651/T;                    %These I don't think need to be reset every time, just recalculated because we aren't addin anything to them
% c_t =  2050; %1925*(T/250); 
% %visc = ;
% 
% 
% depth = [];
% depth(1) = hypot.*sind(sub_A);     %  use sohcahtoa
% Temp(1) = 100;  %
% 
% 
% 
% 
% tau = miu.*rho_i.*g.*depth     ;                 % shear stress
% 
% a = ((2.*tau)./k) ;
% b = (kappa.*depth.*u) ;
% c = (pi.*cosd(phi).*sind(sub_A)) ;
% 
% d_T = a.*((b./c).^(0.5)) ;
% 
% hypot(i) = hypot + 0.04; 
% Temp(i) = Temp + d_T ;         %Something is wrong here maybe?
% 
% 
% disp(Temp)
% 
% 
% 
% i = i+1; 
% 
% end

%% End-Member: No heat diffusion, shell stays at 100 K + friction

% Constants
initial_Temp = 100;             % starting temperature of subducting slab at surface
% k_t = 2.2;%(((0.7).*(651 ./ initial_Temp))+ (0.2 .*0.65) + (0.1.*0));  % + (0.1.*0.65)    % temperature dependent thermal conductivity of ice (Johnson et al 2017)
% c_t = 1925 .* (initial_Temp ./ 250); % temperature dependent thermal capacity of ice (Johnson et al 2017)
g = 1.315;              % gravity of Europa in m/s^2
u = 1.268e-9;           % slab velocity (m/s)  (40 mm/yr)
% C = 2050;               % from temperature of conductive portion of the lid in Kelvin, is linear with depth
% u = 1.268e-9;              % slab velocity (m/s)  (40 mm/yr)
rho_i = 920;            % density of ice in kg/m^3

miu = 0.7 ;            % coefficient of friction (ice sliding on ice)
phi = 0;                % angle between direction of velocity and slab subduction

% Constants for the loop
sub_A = 4.4;             % subduction angle
initial_hypot = 0.04;           % 40 mm, aka moved 40 mm along the hypotenuse to start

 
% Preallocate arrays

max_steps = 5000000;             % number of time steps
depth = zeros(1, max_steps);    % preallocate for depth in meters
Temp = zeros(1, max_steps);     % preallocate for temperature
hypot = zeros(1, max_steps);

hypot(1) = initial_hypot    ;            % intial hypot
depth(1) = hypot(1) .* sind(sub_A);      % initial depth
Temp(1) = initial_Temp;                  % initial temperature

% Adding salt?
addsalt  = 0; 


for z = 1:max_steps-1    % time in years, adjusted to loop from 1 to maxsteps (aka # years I want to simulate)

[k_t] = find_thermalconductivity(Temp(z), 1, 0, 0);
[c_t] = find_heatcapacity(Temp(z), 1, 0, 0);
% [miu] = find_friction(Temp, addsalt);


kappa = k_t ./ (rho_i .* c_t); % thermal diffusivity of water ice
    tau = miu * rho_i * g * depth(z);    % shear stress

    a = ((2 .* tau) ./ k_t);  % assuming k is supposed to be k_t
    b = (kappa .* depth(z) .* u);
    c = (pi .* cosd(phi) .* sind(sub_A));

    d_T = a .* ((b ./ c).^(0.5));
    
    hypot (z + 1) = hypot(z)+ initial_hypot ; 
    depth(z + 1) = hypot(z) .* sind(sub_A);
    Temp(z + 1) = Temp(1) + d_T;
end


plot (Temp, depth./1000, 'LineWidth', 2)


ylabel('Depth (km)')
xlabel('Temperature (K)')
set(gca,'FontSize',15)
set(gca,'ydir','reverse')
% xlim ([0 10]);
% ylim ([0 10]);

%% Endmember: plots of temperature of max case

% Constants
initial_Temp = 100;             % starting temperature of subducting slab at surface
% k_t = 2.2;%(((0.7).*(651 ./ initial_Temp))+ (0.2 .*0.65) + (0.1.*0));  % + (0.1.*0.65)    % temperature dependent thermal conductivity of ice (Johnson et al 2017)
% c_t = 1925 .* (initial_Temp ./ 250); % temperature dependent thermal capacity of ice (Johnson et al 2017)
g = 1.315;              % gravity of Europa in m/s^2
u = 3.171e-9;           % slab velocity (m/s)  (40 mm/yr)
% C = 2050;               % from temperature of conductive portion of the lid in Kelvin, is linear with depth
% u = 1.268e-9;              % slab velocity (m/s)  (40 mm/yr)
rho_i = 920;            % density of ice in kg/m^3


miu = 0.55 ;            % coefficient of friction (ice sliding on ice)
phi = 0;                % angle between direction of velocity and slab subduction

% to find temp of conducting lid with depth 
h_cond = 10000;      % thickness of the conductive portion of the shell (8km here)
T_b = 260;          % Temperature of the base of the conducting lid
depth_intoshell = 1:h_cond ; 

% Constants for the loop

sub_A_values = [4.4, 18.4];  % Example subduction angles
for i = 1:length(sub_A_values)

sub_A = sub_A_values(i);             % subduction angle
initial_hypot = 0.1;           % 40 mm, aka moved 40 mm along the hypotenuse to start

% 
% % Saving files
% saveTemp = 1;
% Temp4degSA = 'Temp6H_Con4degSubA.mat';

%°°°°°°°


% Preallocate arrays

max_steps = 9000000;             % number of time steps
depth = zeros(1, max_steps);    % preallocate for depth in meters
Temp = zeros(1, max_steps);     % preallocate for temperature
hypot = zeros(1, max_steps);

hypot(1) = initial_hypot    ;            % intial hypot
depth(1) = hypot(1) .* sind(sub_A);      % initial depth
Temp(1) =  initial_Temp.*(T_b./initial_Temp).^(depth(1)./h_cond) ;    % initial temperature

% Adding salt?
addsalt  = 0; 


for z = 1:max_steps-1    % time in years, adjusted to loop from 1 to maxsteps (aka # years I want to simulate)


[k_t] = find_thermalconductivity(Temp(z), 0.85, 0.15, 0);
[c_t] = find_heatcapacity(Temp(z), 0.85, 0.15, 0);
% [miu] = find_friction(Temp, addsalt);


kappa = k_t ./ (rho_i .* c_t); % thermal diffusivity of water ice
    tau = miu * rho_i * g * depth(z);    % shear stress

    a = ((2 .* tau) ./ k_t);  % assuming k is supposed to be k_t
    b = (kappa .* depth(z) .* u);
    c = (pi .* cosd(phi) .* sind(sub_A));

    d_T = a .* ((b ./ c).^(0.5));
    
    hypot (z + 1) = hypot(z)+ initial_hypot ; 
    depth(z + 1) = hypot(z) .* sind(sub_A);
    Temp(z + 1) = (initial_Temp.*(T_b./initial_Temp).^(depth(z)./h_cond)) + d_T;   % this is the tricky one, I need to add d_T to whatever the temperature is at that depth
end

% if saveTemp == 1
%     save(Temp4degSA, 'depth', 'Temp') ; 
% end
%

    save(['Temp6H_Con' num2str(sub_A) 'degSubA.mat'], 'depth', 'Temp');
end




hold on


% Plot the reference temperature profile for the conductive lid
T_s = 100;                      % Surface temperature in Kelvin
Temp_shell = T_s .* (T_b ./ T_s) .^ (depth_intoshell ./ h_cond);
plot(Temp_shell, depth_intoshell / 1000, 'LineWidth', 2);

% Loop over each subduction angle
for i = 1:length(sub_A_values)
    % Load the temperature and depth data for each case
    load(['Temp6H_Con' num2str(sub_A_values(i)) 'degSubA.mat'], 'depth', 'Temp');
    
    % Find the index where temperature exceeds 260 K
    stopIndex = find(Temp >= 273, 1, 'first'); % First index where temperature >= 260

    % If stopIndex is found, truncate the data to only include up to that point
    if ~isempty(stopIndex)
        depth = depth(1:stopIndex);     % Truncate depth
        Temp = Temp(1:stopIndex);       % Truncate temperature
    end
    
    % Plot the temperature profile, with depth in km
    plot(Temp, depth / 1000, 'LineWidth', 2);  % Plot each scenario
end

% Add convective portion:

conv_temp = [];
depth_conv_temp = [];

% Finalize the plot
hold off;

xline(204)
xline(273)
ylabel('Depth (km)');
xlabel('Temperature (K)');
set(gca, 'FontSize', 15);
set(gca, 'YDir', 'reverse');  % Reverse the depth axis for better visualization
xlim([100 300]);
ylim([0 15]);


% ylabel('Depth (km)')
% xlabel('Temperature (K)')
% set(gca,'FontSize',15)
% set(gca,'ydir','reverse')
% % xline( 260,'--', 'LineWidth', 2 )
% % fill () 
% xlim ([100 280]);
% ylim ([0 20]);

%% Plotting salts

% Constants
initial_Temp = 100;             % starting temperature of subducting slab at surface
% k_t = 2.2;%(((0.7).*(651 ./ initial_Temp))+ (0.2 .*0.65) + (0.1.*0));  % + (0.1.*0.65)    % temperature dependent thermal conductivity of ice (Johnson et al 2017)
% c_t = 1925 .* (initial_Temp ./ 250); % temperature dependent thermal capacity of ice (Johnson et al 2017)
g = 1.315;              % gravity of Europa in m/s^2
u = 1.268e-9;           % slab velocity (m/s)  (40 mm/yr)
% C = 2050;               % from temperature of conductive portion of the lid in Kelvin, is linear with depth
% u = 1.268e-9;              % slab velocity (m/s)  (40 mm/yr)


miu = 0.55 ;            % coefficient of friction (ice sliding on ice)
phi = 0;                % angle between direction of velocity and slab subduction

% porosity section
phi_i = 0.01;           % initial porosity of subducting slab
g = 1.315;              % gravity of Europa in m/s^2
t_por =  65000000; % time ice shell evolves porosity before subduction (65 Myr)
depth1 = 1:6000;
rho1 = 950;
visc_r = 10^14;         % reference viscosity in Pa*S
T_depth = linspace(100, 260, length(depth1));  % Linearly increasing temperature with depth
Q = 50 ;
T_b = 260 ;
R = 8.314 ; 

P = rho1*g*depth1 ;

visc =  visc_r.*exp((Q./R).* ((1./T_depth)-(1./T_b))) ;

phi_depth = phi_i .* exp((-1).*(t_por.* P)./(visc)) ;

% rho_i = 920;            % density of ice in kg/m^3

% to find temp of conducting lid with depth 
h_cond = 6000;      % thickness of the conductive portion of the shell (8km here)
T_b = 260;          % Temperature of the base of the conducting lid
depth_intoshell = 1:h_cond ; 
% Constants for the loop

sub_A_values = [4.4, 18.4];  % Example subduction angles
for i = 1:length(sub_A_values)

sub_A = sub_A_values(i);             % subduction angle
initial_hypot = 0.04;           % 40 mm, aka moved 40 mm along the hypotenuse to start

% 
% % Saving files
% saveTemp = 1;
% Temp4degSA = 'Temp6H_Con4degSubA.mat';

%°°°°°°°


% Preallocate arrays

max_steps = 3500000;             % number of time steps
depth = zeros(1, max_steps);    % preallocate for depth in meters
Temp = zeros(1, max_steps);     % preallocate for temperature
hypot = zeros(1, max_steps);
rho_i = zeros(1, max_steps);


hypot(1) = initial_hypot    ;            % intial hypot
depth(1) = hypot(1) .* sind(sub_A);      % initial depth
Temp(1) =  initial_Temp.*(T_b./initial_Temp).^(depth(1)./h_cond) ;    % initial temperature
% phi_depth(1) = 
rho_i (1) = 920 - 920.*phi_depth(1) ; 

% Adding salt?
addsalt  = 0; 


for z = 1:max_steps-1    % time in years, adjusted to loop from 1 to maxsteps (aka # years I want to simulate)


[k_t] = find_thermalconductivity(Temp(z), 0.85, 0.15, 0);
[c_t] = find_heatcapacity(Temp(z), 0.85, 0.15, 0);
% [miu] = find_friction(Temp, addsalt);


kappa = k_t ./ (rho_i(z) .* c_t); % thermal diffusivity of water ice
    tau = miu * rho_i(z) * g * depth(z);    % shear stress

    a = ((2 .* tau) ./ k_t); 
    b = (kappa .* depth(z) .* u);
    c = (pi .* cosd(phi) .* sind(sub_A));

    d_T = a .* ((b ./ c).^(0.5));
    
    hypot (z + 1) = hypot(z)+ initial_hypot ; 
    depth(z + 1) = hypot(z) .* sind(sub_A);
    Temp(z + 1) = (initial_Temp.*(T_b./initial_Temp).^(depth(z)./h_cond)) + d_T;   % this is the tricky one, I need to add d_T to whatever the temperature is at that depth

    if z <= length(phi_depth)
    current_phi = phi_depth(z);
else
    current_phi = phi_depth(end);
end
rho_i(z + 1) = 920 - 920 * current_phi;
    
end

% if saveTemp == 1
%     save(Temp4degSA, 'depth', 'Temp') ; 
% end
%

    save(['Temp6H_Con' num2str(sub_A) 'degSubA.mat'], 'depth', 'Temp');
end




hold on


% Plot the reference temperature profile for the conductive lid
T_s = 100;                      % Surface temperature in Kelvin
Temp_shell = T_s .* (T_b ./ T_s) .^ (depth_intoshell ./ h_cond);
plot(Temp_shell, depth_intoshell / 1000, 'LineWidth', 2);

% Loop over each subduction angle
for i = 1:length(sub_A_values)
    % Load the temperature and depth data for each case
    load(['Temp6H_Con' num2str(sub_A_values(i)) 'degSubA.mat'], 'depth', 'Temp');
    
    % Find the index where temperature exceeds 260 K
    stopIndex = find(Temp >= 260, 1, 'first'); % First index where temperature >= 260

    % If stopIndex is found, truncate the data to only include up to that point
    if ~isempty(stopIndex)
        depth = depth(1:stopIndex);     % Truncate depth
        Temp = Temp(1:stopIndex);       % Truncate temperature
    end
    
    % Plot the temperature profile, with depth in km
    plot(Temp, depth / 1000, 'LineWidth', 2);  % Plot each scenario
    
    % % List of target temperatures to display coordinates for
    % target_temps = [204, 236, 252, 260];
    % 
    % % Loop through target temperatures and find the corresponding depth
    % for t = target_temps
    %     % Find the index where temperature exceeds the target value
    %     target_index = find(Temp >= t, 1, 'first');
    % 
    %     if ~isempty(target_index)
    %         % Get the corresponding depth and temperature
    %         target_depth = depth(target_index) / 1000;  % Convert depth to km
    %         target_temp = Temp(target_index);           % Temperature at that depth
    % 
    %         % Plot a marker for the point
    %         plot(target_temp, target_depth, 'ko', 'MarkerFaceColor', 'k');
    % 
    %         % Display the temperature and depth on the plot
    %         text(target_temp + 1, target_depth, sprintf('%.0f, %.2f', target_temp, target_depth), ...
    %             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10, 'Color', 'black');
    %     end
    % end
end

% Add convective portion:
conv_temp = [252 256 257 259 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 261 265 273];
depth_conv_temp = [5.84 5.95 6.0 6.5 7.0 7.5 8 8.5 9.0 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 27.5 28 28.5 29 29.5 30];

plot (conv_temp, depth_conv_temp, 'LineWidth', 2)

xline(204)
xline(236)
xline(252)
xline(260)


% Finalize the plot
hold off;

ylabel('Depth (km)');
xlabel('Temperature (K)');
set(gca, 'FontSize', 15);
set(gca, 'YDir', 'reverse');  % Reverse the depth axis for better visualization
xlim([100 300]);
ylim([0 30]);
addpath('/Users/smenten/Library/CloudStorage/Box-Box/MATLAB/breakyaxis')

breakyaxis([10 25])  % breaks Y-axis between 10 and 25 km






%% Add Porosity into the Model

phi_i = 0.1;           % initial porosity of subducting slab
g = 1.315;              % gravity of Europa in m/s^2
t_por =  65000000; % time ice shell evolves porosity before subduction (65 Myr)
depth1 = 1:6000;
rho1 = 950;
visc_r = 10^14;         % reference viscosity in Pa*S
T_depth = linspace(100, 260, length(depth1));  % Linearly increasing temperature with depth
Q = 50 ;
T_b = 260 ;
R = 8.314 ; 

P = rho1*g*depth1 ;

visc =  visc_r.*exp((Q./R).* ((1./T_depth)-(1./T_b))) ;

phi_depth = phi_i .* exp((-1).*(t_por.* P)./(visc)) ;


plot (phi_depth, depth1./1000, 'LineWidth', 2)

ax = gca;
ylabel('Depth (km)')
xlabel('Porosity')
set(gca,'FontSize',15)
set(gca,'ydir','reverse')
set(gca,'xdir','reverse')
ax.XAxis.Exponent = 0 ; 

% %% Make porosity Temp Dependent (do not need this, above porosity is temp
% dependent
% 
% 
% phi_i = 0.01;            % initial porosity of subducting slab
% g = 1.315;              % gravity of Europa in m/s^2
% t_por =  65000000;      % time ice shell evolves porosity before subduction
% rho1 = 917;
% visc_r = 10^15;         % reference viscosity in Pa*S
% 
% Q = 50 ;
% 
% R = 8.314 ; 
% 
% T_s = 100 ;         % Temperature of surface in Kelvin
% h_cond = 6000;      % thickness of the conductive portion of the shell (6km here)
% T_b = 260;          % Temperature of the base of the conducting lid
% depth_intoshell = 1:h_cond ; 
% 
% sub_A = 14;             % subduction angle
% initial_hypot = 0.04;           % 40 mm, aka moved 40 mm along the hypotenuse to start
% 
% 
% 
% % time dependent porosity, time for for loop!
% 
% 
% % preallocate arrays
% max_steps = 650000;             % number of time steps
% hypot = zeros(1, max_steps);
% depth = zeros(1, max_steps);    % preallocate for depth in meters
% delta_por_time = zeros(1, max_steps);
% P =  zeros(1, max_steps);
% Temp_shell = zeros(1, max_steps);
% visc = zeros(1, max_steps);
% por_time = zeros(1, max_steps);
% 
% 
% hypot(1) = initial_hypot    ;            % intial hypot
% depth(1) = hypot(1) .* sind(sub_A);      % initial depth
% P(1) = rho1*g*depth(1) ;
% Temp_shell(1) = T_s .*(T_b./T_s).^(depth(1)./h_cond) ;
% visc(1) = visc_r.*exp((Q./R).* ((1./Temp_shell(1))-(1./T_b))) ; 
% por_time(1) = phi_i ;
% delta_por_time(1) = (phi_i .* P(1)) /( visc(1));
% % phi_depth = phi_i .* exp((-1).*(t_por.* P)./(visc)) ;
% 
% for z = 1:max_steps-1 
% 
% 
% 
%  % hypot (z + 1) = hypot(z)+ initial_hypot ; 
%  % depth(z + 1) = hypot(z) .* sind(sub_A);
%  % P (z + 1) = rho1.*g.*depth(z) ;
%  % Temp_shell (z + 1) = T_s .*(T_b./T_s).^(depth(z)./h_cond) ;
%  % visc (z + 1) =  visc_r.*exp((Q./R).* ((1./Temp_shell(z))-(1./T_b))) ;
%  % por_time(z + 1) =  ((por_time(z).* P(z))./visc(z)); 
% 
%     hypot(z + 1) = hypot(z) + initial_hypot;
%     depth(z + 1) = hypot(z + 1) .* sind(sub_A);
%     P(z + 1) = rho1 .* g .* depth(z + 1);
%     Temp_shell(z + 1) = T_s .* (T_b ./ T_s) .^ (depth(z + 1) ./ h_cond);
%     visc(z + 1) = visc_r .* exp((Q ./ R) * ((1 ./ Temp_shell(z + 1)) - (1 ./ T_b)));
% 
%     % Update porosity (am i forgetting to compound porosity?)
% 
%       % Update porosity
% delta_por_time(z + 1) = 100.* (por_time(z) .* P(z + 1)) ./ visc(z + 1);
% por_time(z + 1) = por_time(z) - delta_por_time(z + 1);  % Update por_time correctly
% 
% end
% 
% 
% plot (por_time, depth./1000, 'LineWidth', 2)
% 
% ax = gca;
% ylabel('Depth (km)')
% xlabel('Porosity')
% set(gca,'FontSize',15)
% set(gca,'ydir','reverse')
% set(gca,'xdir','reverse')
% ax.XAxis.Exponent = 0 ; 


%% Plot temp of ice shell with depth:

% should be not quite linear with depth because of temp dependent thermal
% conductivity

T_s = 100 ;         % Temperature of surface in Kelvin
h_cond = 8000;      % thickness of the conductive portion of the shell (6km here)
T_b = 260;          % Temperature of the base of the conducting lid
depth_intoshell = 1:h_cond ; 


Temp_shell = T_s .*(T_b./T_s).^(depth_intoshell./h_cond) ;

plot(Temp_shell, depth_intoshell./1000, 'LineWidth', 2)


ylabel('Depth (km)')
xlabel('Temperature (K)')
set(gca,'FontSize',15)
set(gca,'ydir','reverse')

%% Histogram Plot

% Distance = {'0–10', '10–20', '20–30', '30–40', '40–50', '50–60', '60–70', '70–80', '80–90', '90–95'};

Frequency = [2 6 3 4 2 3 1 3 2 1];

bar(Frequency)

set(gca,'xticklabel',{'0–10', '10–20', '20–30', '30–40', '40–50', '50–60', '60–70', '70–80', '80–90', '90–95'});

xlabel('Distance from Subduction Zone (km)')
ylabel('# Mapped Cryolavas')
ylim ([0 10]);



%% Figure 3 Plot

meltdepth10mm_44 = [2.91 4.28 5.58 6.81];
meltdepth10mm_184 = [2.95 4.39 5.78 7.14];

meltdepth40mm_44 = [2.84 4.09 5.22 6.24];
meltdepth40mm_184 = [2.91 4.28 5.6 6.80];

meltdepth100mm_44 = [2.76 3.90 4.89 5.70];
meltdepth100mm_184 = [2.87 4.18 5.38 6.46];

cond_thickness = [4 6 8 10] ;


plot (cond_thickness, meltdepth10mm_44, 'LineWidth', 2, 'color', [0 0.4470 0.7410] ) ;

hold on

plot (cond_thickness, meltdepth10mm_184, 'LineWidth', 2, 'color',  [0 0.4470 0.7410]) ;
plot (cond_thickness, meltdepth40mm_44, 'LineWidth', 2, 'color', [0.8500 0.3250 0.0980] ) ;
plot (cond_thickness, meltdepth40mm_184, 'LineWidth', 2, 'color', [0.8500 0.3250 0.0980]) ;
plot (cond_thickness, meltdepth100mm_44, 'LineWidth', 2, 'color', [0.9290 0.6940 0.1250]) ;
plot (cond_thickness, meltdepth100mm_184, 'LineWidth', 2, 'color', [0.9290 0.6940 0.1250]) ;
hold off

ylabel('Shallowest Depth of Melting (km)')
xlabel('Conductive Lid Thickness (km)')

xlim ([4 10]);
ylim ([2 8]);
set(gca,'ydir','reverse')

%% Plotting salt v2
% Plotting salts

% --- Constants ---
initial_Temp = 100;             % Initial slab surface temperature (K)
g = 1.315;                      % Europa gravity (m/s^2)
u = 1.268e-9;                   % Slab velocity (m/s)
miu = 0.7;                      % Ice-ice friction coefficient
phi = 0;                        % Angle between velocity and slab (deg)

% --- Porosity & depth setup ---
phi_i = 0.01;                   % Initial porosity
t_por = 6.5e7;                  % Time before subduction (s)
depth1 = 1:6000;                % Depth range (m)
rho1 = 950;                     % Ice density (kg/m^3)
T_depth = linspace(100, 260, length(depth1));
Q = 50;                         % Activation energy (J/mol)
T_b = 260;                      % Basal temperature (K)
R = 8.314;                      % Gas constant (J/mol·K)

% --- Viscosity & Porosity profile ---
P = rho1 * g * depth1;
visc_r = 1e14;
visc = visc_r .* exp((Q ./ R) .* ((1 ./ T_depth) - (1 ./ T_b)));
phi_depth_profile = phi_i .* exp((-t_por .* P) ./ visc);

% --- Conductive lid setup ---
h_cond = 6000;                      % Conductive shell thickness (m)
depth_intoshell = 1:h_cond;         % Depth into shell for plotting
T_s = 100;                          % Surface temp (K)
Temp_shell = T_s .* (T_b ./ T_s) .^ (depth_intoshell ./ h_cond);  % Reference temp profile

% --- Subduction angles ---
sub_A_values = [4.4, 18.4];

for i = 1:length(sub_A_values)
    sub_A = sub_A_values(i);
    initial_hypot = 0.04;  % Initial displacement (m)
    max_steps = 3500000;

    % --- Preallocate arrays ---
    depth = zeros(1, max_steps);
    Temp = zeros(1, max_steps);
    hypot = zeros(1, max_steps);
    rho_i_arr = zeros(1, max_steps);

    hypot(1) = initial_hypot;
    depth(1) = hypot(1) * sind(sub_A);
    Temp(1) = initial_Temp .* (T_b / initial_Temp) .^ (depth(1) / h_cond);

    % Interpolate phi_depth for the first point
    phi_z = interp1(depth1, phi_depth_profile, depth(1), 'linear', 'extrap');
    rho_i_arr(1) = 920 - 920 * phi_z;

    for z = 1:max_steps - 1
        % Thermal properties
        [k_t] = find_thermalconductivity(Temp(z), 0.85, 0.15, 0);
        [c_t] = find_heatcapacity(Temp(z), 0.85, 0.15, 0);
        kappa = k_t / (rho_i_arr(z) * c_t);

        % Shear stress and temperature change
        tau = miu * rho_i_arr(z) * g * depth(z);
        a = (2 * tau) / k_t;
        b = kappa * depth(z) * u;
        c = pi * cosd(phi) * sind(sub_A);
        d_T = a * sqrt(b / c);

        % Update values
        hypot(z + 1) = hypot(z) + initial_hypot;
        depth(z + 1) = hypot(z + 1) * sind(sub_A);
        Temp(z + 1) = initial_Temp .* (T_b / initial_Temp) .^ (depth(z + 1) / h_cond) + d_T;

        phi_z = interp1(depth1, phi_depth_profile, depth(z + 1), 'linear', 'extrap');
        rho_i_arr(z + 1) = 920 - 920 * phi_z;

        % Stop if temperature exceeds melting
        if Temp(z + 1) >= 260
            break;
        end
    end

    % Truncate arrays to actual simulation length
    valid_length = z + 1;
    depth = depth(1:valid_length);
    Temp = Temp(1:valid_length);

    % Save data
    save(['Temp6H_Con' num2str(sub_A) 'degSubA.mat'], 'depth', 'Temp');
end

% --- Plotting Section ---

hold on
plot(Temp_shell, depth_intoshell / 1000, 'LineWidth', 2);  % Reference profile

for i = 1:length(sub_A_values)
    load(['Temp6H_Con' num2str(sub_A_values(i)) 'degSubA.mat'], 'depth', 'Temp');
    plot(Temp, depth / 1000, 'LineWidth', 2);
end

% Add convective profile
conv_temp = [252 256 257 259 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 261 265 273];
depth_conv_temp = [5.84 5.95 6.0 6.5 7.0 7.5 8 8.5 9.0 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 27.5 28 28.5 29 29.5 30];
plot(conv_temp, depth_conv_temp, 'LineWidth', 2);

% Highlight key temperature thresholds
xline(204); xline(236); xline(252); xline(260);

% Finalize plot
xlabel('Temperature (K)');
ylabel('Depth (km)');
set(gca, 'YDir', 'reverse', 'FontSize', 15);
xlim([100 300]);
ylim([0 30]);
hold off;

