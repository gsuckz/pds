load('ws.mat'); %Esto para cargar los filtros creados durante el
%desarrollo del programa
%Se cargan los archivos de audio
fs = 8000;
[s16k, fs16] = audioread('16khz.wav');
[s24k, fs24] = audioread('24khz.wav');
[s32k, fs32] = audioread('32khz.wav');

%GRAFICAR SEÑALES (Tiempo y Espectro)
%{
ESPACIO PARA los Graficos
%}

% Adecuación de la señal al estándar de banda base fijado por CCITT/ITU para telefonía fija (300 Hz a 3,4 KHz)
% Filtro Pasa Banda
%Se generan los filtros a través del a herramienta de filtedesign (GUI) 

s1 = filter(f16k1, s16k); 
s2 = filter(f24k1, s24k); 
s3 = filter(f32k1, s32k);


%GRAFICAR SEÑALES (Tiempo y Espectro)
%{
ESPACIO PARA los Graficos
%}

%Se re-samplea (? a 8kHz las señales de banda base

s1_8k = downsample(s1,2);
s2_8k = downsample(s2,3);
s3_8k = downsample(s3,4);

%Falta fijar a 12-bits la resolución

%Se aumenta la cantidad de muestras a 120k
s1_120k1 = upsample(s1_8k,15);
s2_120k1 = upsample(s2_8k,15);
s3_120k1 = upsample(s3_8k,15);

%Igualo la longitud de los vectores

L = max([length(s1_120k1),length(s2_120k1),length(s3_120k1)]);
s1_120k = [s1_120k1 , zeros(1, L - length(s1_120k1))];
s2_120k1 = reshape(s2_120k1, 1, []);
s2_120k = [s2_120k1 , zeros(1, L - length(s2_120k1))];
s3_120k1 = reshape(s3_120k1, 1, []);
s3_120k = [s3_120k1 , zeros(1, L - length(s3_120k1))];









save('ws.mat');