clear; close all; clc;

%% Parámetros
Fs = 2000;        % Frecuencia de muestreo original [Hz]
Ts = 1/Fs;        
N = 2048;         % Número de muestras
t = (0:N-1)*Ts;   % Vector de tiempo

%% Señal de banda base (ejemplo: suma de dos tonos)
f1 = 50; f2 = 120;     % Frecuencias en banda base
x_bb = cos(2*pi*f1*t) + 0.7*sin(2*pi*f2*t);

%% Modulamos en banda lateral superior (traslado a fc)
fc = 400;  
x_usb = x_bb .* exp(1j*2*pi*fc*t);

%% Downsampling
M = 4;                  % Factor de diezmado
y = downsample(x_usb, M);
Fs_ds = Fs/M;           % Nueva frecuencia de muestreo

%% Gráficos de espectro
f = (-N/2:N/2-1)*(Fs/N);       % Eje frecuencia antes
f_ds = (-N/2:N/2-1)*(Fs_ds/N); % Eje frecuencia después

X_usb = fftshift(abs(fft(x_usb, N))/N);
Y = fftshift(abs(fft(y, N))/N);

figure;
subplot(2,1,1);
plot(f, X_usb);
title('Espectro antes del Downsampling (USB)');
xlabel('Frecuencia [Hz]'); ylabel('|X(f)|');

subplot(2,1,2);
plot(f_ds, Y);
title('Espectro después del Downsampling');
xlabel('Frecuencia [Hz]'); ylabel('|Y(f)|');

