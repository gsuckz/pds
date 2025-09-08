clc; clear; close all;

% Archivos de entrada (uno por cada frecuencia de muestreo)
archivos16khz = {'canal1_16_khz.wav', 'canal2_16_khz.wav', 'canal3_16_khz.wav'}; %Incluyo arhivos nombrados por su fsampleo
archivos24khz = {'canal1_24_khz.wav', 'canal2_24_khz.wav', 'canal3_24_khz.wav'}; %Incluyo arhivos nombrados por su fsampleo
archivos32khz = {'canal1_32_khz.wav', 'canal2_32_khz.wav', 'canal3_32_khz.wav'}; %Incluyo arhivos nombrados por su fsampleo
archivos = {archivos16khz,archivos24khz,archivos32khz};
fs_list = [16000, 24000, 32000]; %Genero un vector con las fsamples (puedo hacerlo porque son conocidas)
colores = lines(3);  % Colores para graficar

