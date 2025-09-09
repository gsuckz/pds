
load('archivos_procesados.mat', 'archivos_procesados');
factor = 14;        % Upsampling
%Filtros Multiplexado
fs = 120000;
banda1 = [12300 15400];       % Banda de paso
banda2 = [16300 19400];       % Banda de paso
banda3 = [20300 23400];       % Banda de paso
bandas = {banda1, banda2, banda3};
a = [1 0];               % Magnitudes: 1 en paso, 0 en stop
dev = [0.01 0.01];       % tolerancia: 1% en banda paso y rechazo
filtro = {};

%for k = 1:3
%[N, Fo, Ao, W] = firpmord(banda{k}, a, dev, fs);  % calcula orden y parámetros
%filtro{k} = firpm(N, Fo, Ao, W);                   % diseño FIR óptimo
%end
for k = 1:3
    % Normalizar frecuencias a Nyquist
    Wn = bandas{k}/(fs/2);  % <--- IMPORTANTE: entre 0 y 1
    filtro{k} = fir1(14, Wn, 'bandpass'); %Uso orden 15 para hacer 15 muestras de salida con cada muestra de entrada
end
save('filtro.mat', 'filtro');

canal1 = audioread(archivos_procesados{1});
canal2 = audioread(archivos_procesados{2});
canal3 = audioread(archivos_procesados{3});

L = length(canal1) ; %Todos son iguales, de todas formas esto no sería conocido en un sistema de tiempo real


salidaMux = zeros((L+1)*factor,1);  % inicializar salida
salidaSinFiltrar = zeros((L+1)*factor,1);  % inicializar salida
multDesplazamintoFrecuncia = 1;
for n = 1:L
    % Posición inicial en la señal de salida
    idx = (n-1)*factor + 1;
    muestra1 = canal1(n) * multDesplazamintoFrecuncia;
    muestra2 = canal2(n);
    muestra3 = canal3(n) * multDesplazamintoFrecuncia;
    multDesplazamintoFrecuncia = multDesplazamintoFrecuncia * (-1);
    % Cada ciclo genera 15 muestras
    salidaMux(idx:idx+factor) = muestra1*filtro{1} + muestra2*filtro{2} + muestra3*filtro{3};
    salidaSinFiltrar(idx:idx+factor) = muestra1+ muestra2 + muestra3;

end

save('salidaMux.mat','salidaMux');

%% Señal de salida
fs_salida = 120000;  % frecuencia de muestreo
x = salidaMux;

%% FFT
N = length(x);
X = fft(x);
f = linspace(0, fs_salida, N);

%% Magnitud en dB
magX = 20*log10(abs(X));

%% Graficar espectro completo
figure;
plot(f, magX);
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Espectro de la señal multiplexada a 120 kHz');
grid on;

%% Graficar solo hasta Nyquist (60 kHz)
figure;
plot(f(1:N/2), magX(1:N/2));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Espectro hasta Nyquist (60 kHz)');
grid on;

%% Graficar forma de onda señal multiplexada 
t = (0:N-1)/fs_salida;
figure;
plot(t, x);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Forma de onda de la señal multiplexada a 120 kHz');
grid on;
