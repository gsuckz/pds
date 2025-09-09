%% Parámetros
fs = 120e3;          % Frecuencia de muestreo original
Ts = 1/fs;
N = 2048;            % Cantidad de muestras
t = (0:N-1)*Ts;

%% Señal de prueba: dos senoidales
f1 = 10e3;   % Dentro de la nueva banda de Nyquist
f2 = 35e3;   % Causará aliasing al diezmar
x = cos(2*pi*f1*t) + cos(2*pi*f2*t);

%% FFT helper function
plotFFT = @(sig, Fs, ttl) ...
    plot(linspace(0, Fs, length(sig)), 20*log10(abs(fft(sig))), 'LineWidth',1.2); 

%% Diezmado (M = 2)
M = 2;
y_no_filtro = x(1:M:end);   % Descartar muestras directamente
fs_new = fs/M;              % Nueva frecuencia de muestreo

%% Con filtro antialias
% Pasa bajos con corte en fs_new/2 = 30 kHz
fc = fs_new/2 / (fs/2);    % Normalizado
h = fir1(64, fc, 'low');   % Filtro FIR
x_filtrado = filter(h, 1, x);
y_filtrado = x_filtrado(1:M:end);

%% Graficar FFTs
figure;
subplot(3,1,1);
plotFFT(x, fs, 'Señal original (fs=120 kHz)');
title('Señal original (120 kHz)');
xlabel('Frecuencia (Hz)'); ylabel('Magnitud (dB)'); grid on;

subplot(3,1,2);
plotFFT(y_no_filtro, fs_new, 'Diezmada sin filtro');
title('Diezmada sin filtro (60 kHz)');
xlabel('Frecuencia (Hz)'); ylabel('Magnitud (dB)'); grid on;

subplot(3,1,3);
plotFFT(y_filtrado, fs_new, 'Diezmada con filtro');
title('Diezmada con filtro antialias (60 kHz)');
xlabel('Frecuencia (Hz)'); ylabel('Magnitud (dB)'); grid on;
