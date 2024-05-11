clc

C_m = 1;     
g_Na = 120;   
E_Na = 50;   
g_K = 36;     
E_K = -77;    
g_F = 0.3;    
E_F = -54.4;   


dt = 0.01;                 
T = 50;                     
t = 0:dt:T;                 
pasos = numel(t);

V = zeros(size(t));  
m = zeros(size(t));  
h = zeros(size(t));  
n = zeros(size(t)); 
V(1) = -65;

x=input("Inicio de la corriente: ");
    y=input("Fin de la corriente: ");
    z=input("Cantidad de corriente aplicada (uA): ");

    I = zeros(size(t));  
    I(find(t >= x & t <= y)) = z;  

for i = 1:pasos-1

alpha_m = (0.1*(V(i)+40)) / (1 - exp(-(V(i)+40)/10));
beta_m = 4 * exp(-(V(i)+65)/18);
alpha_h = 0.07 * exp(-(V(i)+65)/20);
beta_h = 1 / (1 + exp(-(V(i)+35)/10));
alpha_n = (0.01*(V(i)+55)) / (1 - exp(-(V(i)+55)/10));
beta_n = 0.125 * exp(-(V(i)+65)/80);
    

    iNa = g_Na * m(i)^3 * h(i) * (V(i) - E_Na);
    iK = g_K * n(i)^4 * (V(i) - E_K);
    iF = g_F * (V(i) - E_F);
     
    dvdt = (1/C_m) * (I(i)-iNa - iK - iF);
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

