% CIV102 Bridge Design Project

clear; close all;

% Step 0: Initialize Parameters
L = 1200;   % Length of the bridge
n = 1200;   % Discretize into 1 mm segments
P = 400;    % Total weight of the train [N]
x = linspace(0, L, n+1);    % x-axis



% Step 1: SFD, BMD Under Train Loading
x_train = [52, 228, 392, 568, 732, 908]; % train load locations
P_train = [1 1 1 1 1 1] * (P/6);

n_train = 25;    % number of train start locations
SFDi = zeros(n_train, n+1);     % 1 SFD for each train start location
BMDi = zeros(n_train, n+1);     % 1 BMD for each train start location

for i = 1:n_train
    x_start_location = 10 * i;  % start location of the train, acts like an adjustment factor
    x_train_adjusted = x_train + x_start_location;  % adjusted train load locations

    % sum of moments at A equation, rearranged for R_y
    R_y = (sum(P_train .* x_train_adjusted)) / L;

    % sum of F_y equation, rearranged for L_y
    L_y = sum(P_train) - R_y;


    v = L_y; % initial shear force at left support
    SFDi(i, 1:x_train_adjusted(1)) = v; % filling in a shear force value for a given range of x-values
    v = v - P_train(1); % updating the shear force
    SFDi(i, x_train_adjusted(1):x_train_adjusted(2)) = v; % filling in the updated value for a given range of x-values
    v = v - P_train(2);
    SFDi(i, x_train_adjusted(2):x_train_adjusted(3)) = v;
    v = v - P_train(3);
    SFDi(i, x_train_adjusted(3):x_train_adjusted(4)) = v;
    v = v - P_train(4);
    SFDi(i, x_train_adjusted(4):x_train_adjusted(5)) = v;
    v = v - P_train(5);
    SFDi(i, x_train_adjusted(5):x_train_adjusted(6)) = v;
    v = v - P_train(6);
    SFDi(i, x_train_adjusted(6):n) = v;
    v = v + R_y;    % final shear force update at the right support
    SFDi(i, n+1) = v;
    BMDi(i,:) = cumtrapz(SFDi(i,:));   % integrating SFDi
end

SFD = max(abs(SFDi));   % the absolute maximum values of shear force at all locations
BMD = max(BMDi);        % the maximum value of all bending moments at all locations
disp(BMDi(n+1))


% Plotting the SFD and BMD
figure(1)
hold on
%plot(x, SFDi(1,:), "k")
plot(x, SFDi, 'k')  % SFD for individual location tests
plot(x, SFD, "r")   % SFD for maximum shear force
title("Shear Force Diagrams")
xlabel('Distance Along Bridge (mm)')
ylabel('Shear Force (N)')
hold off

figure(2)
hold on
%plot(x, BMDi(1,:), "k")
plot(x, BMDi, "k")  % BMD for individual location tests
plot(x, BMD, "r")   % BMD for maximum bending moment
title("Bending Moment Diagrams")
xlabel('Distance Along Bridge (mm)')
ylabel('Bending Moment (N m)')
hold off
