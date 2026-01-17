clc; clear; close all; 

% Time array from the image
time_arr = { ...
'10:00';'10:10';'10:20';'10:30';'10:40';'10:50';...
'11:00';'11:10';'11:20';'11:30';'11:40';'11:50';...
'12:00';'12:10';'12:20';'12:30';'12:40';'12:50';...
'13:00';'13:10';'13:20';'13:30';'13:40';'13:50';...
'14:00';'14:10';'14:20';'14:30';'14:40';'14:50';...
'15:00';'15:10';'15:20';'15:30';'15:40';'15:50';...
'16:00'};

% Time calculations
dt = 10; % minutes interval
time_tr_min = (0:dt:(dt*(numel(time_arr)-1)))';
time_tr_s = time_tr_min * 60;

% Water Temperature (Tw)
Tw_transient = [ 
26.0; 29.1; 32.5; 36.0; 39.6; 43.2; 46.8; 50.5; 54.2; 57.8;...
61.4; 65.0; 68.6; 72.1; 75.7; 79.1; 82.3; 85.2; 87.9;...
90.4; 92.6; 95.0; 95.0; 95.0; 95.0; 95.0; 94.1; 94.1; 94.1;...
93.3; 93.3; 92.5; 92.5; 91.7; 91.7; 90.9; 89.3];

% Irradiance (I)
I_tr = [
555; 574; 593; 617; 629; 641; 648; 659; 676; 687;...
703; 710; 720; 725; 750; 764; 778; 778; 778;...
762; 746; 727; 727; 711; 708; 700; 687; 665; 655;...
637; 616; 597; 585; 572; 548; 528; 495];

% Ambient Temperature (Ta)
Ta_tr = [
26.0; 26.0; 27.0; 27.0; 28.0; 29.0; 29.0; 30.0; 29.5; 30.5;...
30.0; 30.5; 32.0; 32.0; 32.0; 31.5; 31.0; 32.0; 33.0;...
32.5; 33.0; 32.0; 32.0; 32.0; 31.5; 33.0; 32.5; 33.0; 33.0;...
32.5; 33.0; 34.0; 33.5; 33.5; 34.0; 33.0; 32.5];

%% ----- Constants -----
M  = 2; Cp = 4186; a = 0.48; b = 0.42; h = 0.5;

pI = polyfit(time_tr_s, I_tr, 2);
pTa = polyfit(time_tr_s, Ta_tr, 2);

t_fit = linspace(min(time_tr_s), max(time_tr_s), 200);
I_fit = polyval(pI, t_fit);
Ta_fit = polyval(pTa, t_fit);
x1 = pI(1); x2 = pI(2); x3 = pI(3);
y1 = pTa(1); y2 = pTa(2); y3 = pTa(3);
T0 = 26;

I = @(t) x1.*(t.^2) + x2.*t + x3;
Ta = @(t) y1.*(t.^2) + y2*t + y3;
%% ----- Phi Loop -----
phi_range = 40:1:240;
F2_avg_arr = zeros(size(phi_range));

for k = 1:length(phi_range)
    phi_deg = phi_range(k);
    
    hi = h .* sqrt(2 .* (1 - cosd(phi_deg)));
    Ai = hi * (a + b) / 2;
     
    tau_alpha = 0.4752;
    UL = 4.8287;
    F1_pred = tau_alpha / UL;
    

    [t, Tw_calc] = ode45(@(t, T) proj(t, Ai, tau_alpha, UL, M, Cp, T, x1, x2, x3, y1, y2, y3), time_tr_s, T0);

    % ----- Compute F2 (Predicted Tw) -----
    % Tw_calc = interp1(t, Tw_calc, time_tr_s, 'spline');
    n = length(Tw_calc);
    F2_pred = NaN(n-1,1);

    for i = 1:n-1
        t1 = t(i);  t2 = t(i+1);
        Tw1 = Tw_calc(i);  Tw2 = Tw_calc(i+1);
        Ta1 = Ta(t1);    Ta2 = Ta(t2);
        Ta_avg = 0.5*(Ta2+Ta1);
        Is = 0.5*(I(t2)+I(t1));
        dt = t2-t1;
        num = 1 - (1/F1_pred)*((Tw1 - Ta_avg)/Is);
        den = 1 - (1/F1_pred)*((Tw2 - Ta_avg)/Is);
        if num > 0 && den > 0
            F2_pred(i) = F1_pred*(M*Cp)/(Ai*dt) * log(num/den);
        end
    end

    % ----- Mean F2 via Integration -----
    valid_pred = ~isnan(F2_pred);
    t_pred = time_tr_s(1:end-1);
    t_pred = t_pred(valid_pred);
    F2_pred_valid = F2_pred(valid_pred);
    mean_F2_pred_int = trapz(t_pred, F2_pred_valid) / trapz(t_pred, ones(size(F2_pred_valid)));

    F2_avg_arr(k) = mean_F2_pred_int;
end

%% ----- Results -----
[max_F2, idx_max] = max(F2_avg_arr);
phi_opt = phi_range(idx_max);

fprintf('\n===== RESULTS OVER PHI RANGE =====\n');
fprintf('Maximum avg F2 = %.6f at phi = %.2f°\n', max_F2, phi_opt);

F2_avg_smooth = smoothdata(F2_avg_arr, 'movmean', 5);
% Plot F2 vs Phi
figure;
plot(phi_range, F2_avg_arr,LineWidth=2);
xlabel('\phi (deg)');
ylabel('Average F2');
title('Variation of Mean F2 with Reflector Angle \phi');
grid on;


function dTdt = proj(t, Ai, tau_alpha, UL, M, Cp, T,x1, x2, x3, y1, y2, y3)
    I = x1*t^2+x2*t+x3;
    Ta = y1*t^2+y2*t+y3;
    A1 = M*Cp/Ai;
    A2 = tau_alpha;
    A3 = UL;   
    dTdt = (A2*I-A3*(T-Ta))/A1;
end
