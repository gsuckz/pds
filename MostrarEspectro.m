clc; clear; close all;

% Archivos de entrada (uno por cada frecuencia de muestreo)
archivos16khz = {'canal1_16khz.wav', 'canal2_16khz.wav', 'canal3_16khz.wav'}; %Incluyo arhivos nombrados por su fsampleo
archivos24khz = {'canal1_24_khz.wav', 'canal2_24_khz.wav', 'canal3_24_khz.wav'}; %Incluyo arhivos nombrados por su fsampleo
archivos32khz = {'canal1_32_khz.wav', 'canal2_32_khz.wav', 'canal3_32_khz.wav'}; %Incluyo arhivos nombrados por su fsampleo
archivos = {archivos16khz,archivos24khz,archivos32khz};
fs_list = [16000, 24000, 32000]; %Genero un vector con las fsamples (puedo hacerlo porque son conocidas)
colores = lines(3);  % Colores para graficar

%% ADQUISICI�N Y ESPECTRO
fprintf('--- An�lisis espectral de se�ales originales ---\n')
for i= 1:3
for k = 1:3 %Bucle para leer las 3 se�ales
    [x, fs] = audioread(archivos{i}{k});  %Leo la se�al obetniendo el vector y su fsamle
    t = (0:length(x)-1)/fs; %genero el Vector tiempo
    % Espectro
    N = length(x);
    X = abs(fft(x));
    f = linspace(0, fs, N); %Genero el espectro hasta Fsample 
    figure(i);
    subplot(3,1,k);
    plot(f(1:N/2), 20*log10(X(1:N/2)), 'Color', colores(k,:));
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title([archivos{i}{k} ' - fs = ' num2str(fs/1000) ' kHz']);
    grid on;
end
end