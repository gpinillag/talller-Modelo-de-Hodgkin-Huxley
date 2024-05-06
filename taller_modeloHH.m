clc
% Parámetros del modelo de Hodgkin-Huxley
C_m = 1;      % Capacitancia (uF/cm^2)
g_Na = 120;   % Conductancia de sodio (mS/cm^2)
E_Na = 115;   % Potencial de equilibrio del sodio (mV)
g_K = 36;     % Conductancia de potasio (mS/cm^2)
E_K = -12;    % Potencial de equilibrio del potasio (mV)
g_L = 0.3;    % Conductancia de fuga (mS/cm^2)
E_L = 10.6;   % Potencial de equilibrio de la fuga (mV)

dt = 0.01;                  % Paso de tiempo (ms)
T = 50;                     % Duración de la simulación (ms)
t = 0:dt:T;                 % Vector de tiempo
num_steps = numel(t);
V = zeros(1, num_steps);    % Potencial de membrana (mV)
m = zeros(1, num_steps);    % Variable de activación del sodio
n = zeros(1, num_steps);    % Variable de activación del potasio
h = zeros(1, num_steps);    % Variable de inactivación del sodio

% Condiciones iniciales
V(1) = -65;     % Potencial de membrana inicial (mV)
m(1) = 0.5;     % Valor inicial de m
h(1) = 0.06;    % Valor inicial de h
n(1) = 0.5;     % Valor inicial de n

% Simulación del modelo de Hodgkin-Huxley
for i = 1:num_steps-1

    alpha_m = 0.1 * (V(i) + 40) / (1 - exp(-(V(i) + 40) / 10));
    beta_m = 4 * exp(-(V(i) + 65) / 18);
    alpha_h = 0.07 * exp(-(V(i) + 65) / 20);
    beta_h = 1 / (1 + exp(-(V(i) + 35) / 10));
    alpha_n = 0.01 * (V(i) + 55) / (1 - exp(-(V(i) + 55) / 10));
    beta_n = 0.125 * exp(-(V(i) + 65) / 80);
    
    % Ecuaciones diferenciales del modelo de Hodgkin-Huxley
    dvdt = (1/C_m) * (-g_Na * m(i)^3 * h(i) * (V(i) - E_Na) - g_K * n(i)^4 * (V(i) - E_K) - g_L * (V(i) - E_L));
    dmdt = alpha_m * (1 - m(i)) - beta_m * m(i);
    dndt = alpha_n * (1 - n(i)) - beta_n * n(i);
    dhdt = alpha_h * (1 - h(i)) - beta_h * h(i);
    
    % Actualización de las variables
    V(i+1) = V(i) + dt * dvdt;
    m(i+1) = m(i) + dt * dmdt;
    n(i+1) = n(i) + dt * dndt;
    h(i+1) = h(i) + dt * dhdt;
end

% Graficar el potencial de acción a lo largo del tiempo
figure;
plot(t, V);
xlabel('Tiempo (ms)');
ylabel('Potencial de membrana (mV)');
title('Modelo de Hodgkin-Huxley PA');
grid on;

