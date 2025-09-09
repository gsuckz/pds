%% Señal de ejemplo multiplexada
fs = 120000;           % Frecuencia de muestreo
t = 0:1/fs:0.01;       % 10 ms de señal

% Dos señales en banda base (fs = 8kHz antes de interpolar)
s1 = sin(2*pi*500*t(1:15:end));   % Señal 1 (500 Hz en 8 kHz)
s2 = cos(2*pi*1000*t(1:15:end));  % Señal 2 (1000 Hz en 8 kHz)

% Interpolamos manualmente (upsampling factor 15)
s1_up = upsample(s1,15);
s2_up = upsample(s2,15);

% Filtros pasabanda para multiplexar
b1 = fir1(64,[0.1 0.2]);   % banda baja (~12–24 kHz, simplificado)
b2 = fir1(64,[0.3 0.4]);   % banda más alta (~36–48 kHz)

x1 = filter(b1,1,s1_up);
x2 = filter(b2,1,s2_up);

% Señal multiplexada
salidaMux = x1 + x2;

%% -------- DEMULTIPLEXADO --------
factor = 15;   % decimación
N = length(salidaMux);

% Filtros pasabanda para recuperar
h1 = b1;   % mismo que en TX
h2 = b2;

% Buffers de salida
canal1 = [];
canal2 = [];

% Procesamiento muestra a muestra
for n = 1:factor:N-factor
    % Tomar 15 muestras de la señal multiplexada
    bloque = salidaMux(n:n+factor-1);
    
    % Aplicar filtros pasa banda
    y1 = sum(bloque .* fliplr(h1(1:factor)));  % salida instantánea canal1
    y2 = sum(bloque .* fliplr(h2(1:factor)));  % salida instantánea canal2
    
    % Guardar la muestra ya decimada
    canal1(end+1) = y1;
    canal2(end+1) = y2;
end

%% Visualización
fs_out = fs/factor;  % 8 kHz
t_out = (0:length(canal1)-1)/fs_out;

figure;
subplot(2,1,1);
plot(t_out, canal1);
title('Canal 1 demultiplexado (8 kHz)');
xlabel('Tiempo (s)'); ylabel('Amplitud');

subplot(2,1,2);
plot(t_out, canal2);
title('Canal 2 demultiplexado (8 kHz)');
xlabel('Tiempo (s)'); ylabel('Amplitud');
