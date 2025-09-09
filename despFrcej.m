%% Parámetros
fs = 120e3;          % Frecuencia de muestreo (Hz)
N = 1024;            % Número de muestras
f0 = 10e3;           % Frecuencia de la señal original (Hz)
n = 0:N-1;           % Índice temporal

%% Señal original: una senoide en 10 kHz
x = cos(2*pi*f0*n/fs);

%% Señal modulada por (-1)^n
x_mod = x .* ((-1).^n);

%% FFT de ambas señales
X = fftshift(fft(x));
X_mod = fftshift(fft(x_mod));

f = linspace(-fs/2, fs/2, N);  % eje de frecuencias

%% Graficar espectros
figure;
subplot(2,1,1);
plot(f, abs(X));
title('Espectro de la señal original (10 kHz)');
xlabel('Frecuencia (Hz)');
ylabel('|X(f)|');
grid on;

subplot(2,1,2);
plot(f, abs(X_mod));
title('Espectro después de multiplicar por (-1)^n');
xlabel('Frecuencia (Hz)');
ylabel('|X_{mod}(f)|');
grid on;
