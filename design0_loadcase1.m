% CIV102 Bridge Design Project

clear; close all;

% Step 0: Initialize Parameters
L = 1200;   % Length of the bridge
n = 1200;   % Discretize into 1 mm segments
P = 400;    % Total weight of the train [N]
x = linspace(0, L, n+1);    % x-axis

% Step 1: SFD, BMD Under Train Loading
load_locations = [0 176 340 516 680 856];    % train load locations
point_loads = [1 1 1 1 1 1] * (P/6);

num_locations = L - load_locations(6) + 1;    % number of train start locations
SFDi = zeros(num_locations, n+1);     % 1 SFD for each train start location
BMDi = zeros(num_locations, n+1);     % 1 BMD for each train start location

for i = 1:num_locations
    adjusted_load_locations = load_locations + i;  % adjusted train load locations
    % sum of moments at A equation, rearranged for R_y (right support force)
    R_y = (sum(point_loads .* (adjusted_load_locations - 1))) / L;

    % sum of F_y equation, rearranged for L_y (left support force)
    L_y = sum(point_loads) - R_y;

    V = L_y; % initial shear force at left support
    SFDi(i, 1:adjusted_load_locations(1)) = V; % filling in a shear force value for a given range of x-values
    V = V - point_loads(1); % updating the shear force
    SFDi(i, (adjusted_load_locations(1)):adjusted_load_locations(2)) = V; % filling in the updated value for a given range of x-values
    V = V - point_loads(2);
    SFDi(i, adjusted_load_locations(2):adjusted_load_locations(3)) = V;
    V = V - point_loads(3);
    SFDi(i, adjusted_load_locations(3):adjusted_load_locations(4)) = V;
    V = V - point_loads(4);
    SFDi(i, adjusted_load_locations(4):adjusted_load_locations(5)) = V;
    V = V - point_loads(5);
    SFDi(i, adjusted_load_locations(5):adjusted_load_locations(6)) = V;
    V = V - point_loads(6);
    SFDi(i, adjusted_load_locations(6):n) = V;
    V = V + R_y;    % final shear force update
    SFDi(i, n+1) = V;
    SFDi = round(SFDi, 4);
    BMDi(i,:) = cumsum(SFDi(i,:));   % integrating SFDi
    BMDi = round(BMDi, 4);
end

SFE = max(abs(SFDi));   % the absolute maximum values of shear force at all locations
BME = max(BMDi);        % the maximum value of all bending moments at all locations

%% 2. Define Bridge Parameters



%{
% = xc, bft, tft,
param = [0, 100, 1.27,...
400, 100, 1.27,...
800, 100, 1.27,...
L, 100, 1.27, ...]
%x_c Location, x, of cross-section change
%bft Top Flange Width
%tft Top Flange Thickness
% Extracting user input assuming linear relationship
bft = interp1(param(:,1), param(:,2), x);
tft = interp1(param(:,1), param(:,3), x);
%% 3. Calculate Sectional Properties
% ybar. location of centroidal axis from the bottom
ybar =
ybot =
ytop =
% I
I =
% Q at centroidal axes
Qcent =
% Q at glue location
Qglue =
%% 4. Calculate Applied Stress
S_top =
S_bot =
T_cent =
T_glue =
%% 5. Material and Thin Plate Buckling Capacities
E = 4000;
mu = 0.2;
S_tens =
S_comp =
T_max =
T_gmax =
S_buck1 =
S_buck2 =
S_buck3 =
T_buck =
%% 6. FOS
FOS_tens =
FOS_comp =
FOS_shear =
FOS_glue =
FOS_buck1 =
FOS_buck2 =
FOS_buck3 =
FOS_buckV =
%% 7. Min FOS and the failure load Pfail
minFOS =
Pf =
%% 8. Vfail and Mfail
Mf_tens =
Mf_comp =
Vf_shear =
Vf_glue =
Mf_buck1 =
Mf_buck2 =
Mf_buck3 =
Vf_buckV =
%}

%% 9. Output plots of Vfail and Mfail
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

subplot(2,3,1) % Shear Failure
hold on; grid on; grid minor;
title("Shear Force Diagram vs. Shear Force Capacities")
%plot(x, Vf_shear, 'r')
%plot(x, -Vf_shear.* SFD, 'r')
plot(x, SFDi(216,:), 'k');
plot(x, SFE, '-r');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
legend('Matboard Shear Failure')
xlabel('Distance along bridge (mm)')
ylabel('Shear Force (N)')


subplot(2,3,2) % Shear Failure
hold on; grid on; grid minor;
plot(x, BME, 'k');
[~, max_moment_idx] = max(BME);  % find index of max moment
plot(max_moment_idx, BME(max_moment_idx), "bo") % plot max moment
plot(x, BMDi(216, :), '-r')
xlabel('Distance along bridge (mm)')




























%{
% Plotting the SFD and BMD
figure(1)
hold on
plot(x, SFDi(2,:), "k")
%plot(x, SFDi, 'k')  % SFD for individual location tests
plot(x, SFD, "r")   % SFD for maximum shear force
title("Shear Force Diagrams")
xlabel('Distance Along Bridge (mm)')
ylabel('Shear Force (N)')
hold off

figure(2)
hold on
plot(x, BMDi(1,:), "k")
%plot(x, BMDi, "k")  % BMD for individual location tests
plot(x, BMD, "r")   % BMD for maximum bending moment
title("Bending Moment Diagrams")
xlabel('Distance Along Bridge (mm)')
ylabel('Bending Moment (N mm)')
[~, max_moment_idx] = max(BMD);  % find index of max moment
plot(max_moment_idx, BMD(max_moment_idx), "bo") % plot max moment
hold off

%}
