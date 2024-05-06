clc

C_m = 1;     
g_Na = 120;   
E_Na = 115;   
g_K = 36;     
E_K = -12;    
g_L = 0.3;    
E_L = 10.6;   

dt = 0.01;                 
T = 50;                     
t = 0:dt:T;                 
num_steps = numel(t);
V = zeros(1, num_steps);    
m = zeros(1, num_steps);    
n = zeros(1, num_steps);    
h = zeros(1, num_steps);    

% Condiciones iniciales
V(1) = -65;     
m(1) = 0.5;     
h(1) = 0.06;    
n(1) = 0.5;     


for i = 1:num_steps-1

    alpha_m = 0.1 * (V(i) + 40) / (1 - exp(-(V(i) + 40) / 10));
    beta_m = 4 * exp(-(V(i) + 65) / 18);
    alpha_h = 0.07 * exp(-(V(i) + 65) / 20);
    beta_h = 1 / (1 + exp(-(V(i) + 35) / 10));
    alpha_n = 0.01 * (V(i) + 55) / (1 - exp(-(V(i) + 55) / 10));
    beta_n = 0.125 * exp(-(V(i) + 65) / 80);
    
    
    dvdt = (1/C_m) * (-g_Na * m(i)^3 * h(i) * (V(i) - E_Na) - g_K * n(i)^4 * (V(i) - E_K) - g_L * (V(i) - E_L));
    dmdt = alpha_m * (1 - m(i)) - beta_m * m(i);
    dndt = alpha_n * (1 - n(i)) - beta_n * n(i);
    dhdt = alpha_h * (1 - h(i)) - beta_h * h(i);
    
    
    V(i+1) = V(i) + dt * dvdt;
    m(i+1) = m(i) + dt * dmdt;
    n(i+1) = n(i) + dt * dndt;
    h(i+1) = h(i) + dt * dhdt;
end


figure;
plot(t, V);
xlabel('Tiempo (ms)');
ylabel('Potencial de membrana (mV)');
title('Modelo de Hodgkin-Huxley PA');
grid on;



