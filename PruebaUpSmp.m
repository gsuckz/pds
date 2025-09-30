clear; clc; close all;

%% Parámetros
fs_in = 8000;        % Frecuencia de muestreo inicial
fs_out = 120000;     % Frecuencia de muestreo deseada
L = fs_out / fs_in;  % Factor de upsampling (15)
t_in = 0:1/fs_in:0.01;   % Tiempo para la señal original (10 ms aprox)

%% Señal de entrada (dos tonos)
f1 = 500;   % Frecuencia del primer tono
f2 = 1200;  % Frecuencia del segundo tono
x = sin(2*pi*f1*t_in) + 0.5*sin(2*pi*f2*t_in);

%% Upsampling (insertar 14 ceros)
x_up = upsample(x, L);   % Inserta 14 ceros entre muestras
t_out = (0:length(x_up)-1)/fs_out;

%% Espectros
Nfft = 4096;
f_in = (-Nfft/2:Nfft/2-1)*(fs_in/Nfft);
X_in = fftshift(abs(fft(x, Nfft))/length(x));

f_out = (-Nfft/2:Nfft/2-1)*(fs_out/Nfft);
X_out = fftshift(abs(fft(x_up, Nfft))/length(x_up));

%% Gráficas
figure('Name','Upsampling Ejemplo','NumberTitle','off');

subplot(2,2,1);
plot(t_in*1000, x, 'b','LineWidth',1.2);
xlabel('Tiempo [ms]'); ylabel('Amplitud');
title('Señal original (fs = 8 kHz)');
grid on;

subplot(2,2,2);
plot(t_out*1000, x_up, 'r','LineWidth',1.2);
xlabel('Tiempo [ms]'); ylabel('Amplitud');
title('Señal upsampled (fs = 120 kHz)');
xlim([0 10]); % primeros 10 ms
grid on;

subplot(2,2,3);
plot(f_in, X_in,'b','LineWidth',1.2);
xlabel('Frecuencia [Hz]'); ylabel('|X(f)|');
title('Espectro señal original');
xlim([-4000 4000]); grid on;

subplot(2,2,4);
plot(f_out, X_out,'r','LineWidth',1.2);
xlabel('Frecuencia [Hz]'); ylabel('|X(f)|');
title('Espectro señal upsampled');
xlim([-60000 60000]); grid on;
