clc; clear; close all;

% --- Time array from the image ---
time_arr = { ...
'10:00';'10:10';'10:20';'10:30';'10:40';'10:50';...
'11:00';'11:10';'11:20';'11:30';'11:40';'11:50';...
'12:00';'12:10';'12:20';'12:30';'12:40';'12:50';...
'13:00';'13:10';'13:20';'13:30';'13:40';'13:50';...
'14:00';'14:10';'14:20';'14:30';'14:40';'14:50';...
'15:00';'15:10';'15:20';'15:30';'15:40';'15:50';...
'16:00'};

% --- Time calculations ---
dt = 10; % minutes interval
time_tr_min = (0:dt:(dt*(numel(time_arr)-1)))';   % minutes
time_tr_hr  = 10 + (time_tr_min / 60);            % start at hour 10
time_tr_s   = time_tr_min * 60;                   % seconds

% --- Experimental Data ---
Tw_transient = [ 
26.0; 29.1; 32.5; 36.0; 39.6; 43.2; 46.8; 50.5; 54.2; 57.8;...
61.4; 65.0; 68.6; 72.1; 75.7; 79.1; 82.3; 85.2; 87.9;...
90.4; 92.6; 95.0; 95.0; 95.0; 95.0; 95.0; 94.1; 94.1; 94.1;...
93.3; 93.3; 92.5; 92.5; 91.7; 91.7; 90.9; 89.3];

I_tr = [
555; 574; 593; 617; 629; 641; 648; 659; 676; 687;...
703; 710; 720; 725; 750; 764; 778; 778; 778;...
762; 746; 727; 727; 711; 708; 700; 687; 665; 655;...
637; 616; 597; 585; 572; 548; 528; 495];

Ta_tr = [
26.0; 26.0; 27.0; 27.0; 28.0; 29.0; 29.0; 30.0; 29.5; 30.5;...
30.0; 30.5; 32.0; 32.0; 32.0; 31.5; 31.0; 32.0; 33.0;...
32.5; 33.0; 32.0; 32.0; 32.0; 31.5; 33.0; 32.5; 33.0; 33.0;...
32.5; 33.0; 34.0; 33.5; 33.5; 34.0; 33.0; 32.5];

%% ----- Constants -----
M  = 2; 
Cp = 4186; 
a = 0.48; 
b = 0.42; 
h = 0.5;

tau_alpha = 0.4752;
UL = 4.8287;
F1_pred = tau_alpha / UL;
fprintf('Initial tau_alpha = %.6g, UL = %.6g, F1 = %.6g\n', tau_alpha, UL, F1_pred);

% Polynomial fit for I and Ta
pI = polyfit(time_tr_s, I_tr, 2);
pTa = polyfit(time_tr_s, Ta_tr, 2);

x1 = pI(1); x2 = pI(2); x3 = pI(3);
y1 = pTa(1); y2 = pTa(2); y3 = pTa(3);

T0 = 26;  % Initial temperature

% --- Tilt angles for simulation ---
phi_values = [30, 40, 60, 75, 90, 105, 120, 180];  % degrees

figure;
plot(time_tr_hr, Tw_transient, 'ko', 'MarkerFaceColor', 'y', 'DisplayName', 'Experimental'); 
hold on;

for i = 1:length(phi_values)
    phi_deg = phi_values(i);
    hi = h * sqrt(2 * (1 - cosd(phi_deg)));
    Ai = hi * (a + b) / 2;

    [t1, T_model] = ode45(@(t, T) proj(t, Ai, tau_alpha, UL, M, Cp, T, x1, x2, x3, y1, y2, y3), time_tr_s, T0);

    % Convert ODE time to actual clock hour (10 → 16)
    t1_hr = 10 + (t1 / 3600);

    label_str = sprintf('\\phi = %d° (h_i = %.3f m)', phi_deg, hi);

    plot(t1_hr, T_model, '-', 'LineWidth', 2, 'DisplayName', label_str);
end

xlabel('Time (Hour of Day)');
ylabel('Water Temperature (°C)');
legend('Location', 'best');
title('Effect of Tilt Angle (\phi) on Water Temperature');
grid on;

% ---- ODE ----
function dTdt = proj(t, Ai, tau_alpha, UL, M, Cp, T, x1, x2, x3, y1, y2, y3)
    I = x1*t^2 + x2*t + x3;
    Ta = y1*t^2 + y2*t + y3;
    A1 = M * Cp / Ai;
    dTdt = (tau_alpha * I - UL * (T - Ta)) / A1;
end
