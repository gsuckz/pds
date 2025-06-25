clc; clear; close all;

% Archivos de entrada (uno por cada frecuencia de muestreo)
archivos = {'16khz.wav', '24khz.wav', '32khz.wav'}; %Incluyo arhivos nombrados por su fsampleo
fs_list = [16000, 24000, 32000]; %Genero un vector con las fsamples (puedo hacerlo porque son conocidas)
colores = lines(3);  % Colores para graficar

%% ADQUISICIÓN Y ESPECTRO
fprintf('--- Análisis espectral de señales originales ---\n')
for k = 1:3 %Bucle para leer las 3 señales
    [x, fs] = audioread(archivos{k});  %Leo la señal obetniendo el vector y su fsamle
    t = (0:length(x)-1)/fs; %genero el Vector tiempo

    % Espectro
    N = length(x);
    X = abs(fft(x));
    f = linspace(0, fs, N); %Genero el espectro hasta Fsample 

    figure(1);
    subplot(3,1,k);
    plot(f(1:N/2), 20*log10(X(1:N/2)), 'Color', colores(k,:));
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title([archivos{k} ' - fs = ' num2str(fs/1000) ' kHz']);
    grid on;
end

%% FILTRO ANTIALIASING ANALÓGICO (para fs = 8 kHz)
%% Diseño del filtro analogico necesario (No implementado, solo se grafica su función de Transferencia)
fp = 3400; Ap = 1;
fs_alias = 4000; As = 40;

fprintf('\n--- Filtro antialiasing analógico ---\n')
[na, Wn_a] = buttord(2*pi*fp, 2*pi*fs_alias, Ap, As, 's');
[ba, aa] = butter(na, Wn_a, 's');
Ha = tf(ba, aa);

figure;
bode(Ha);
title('Filtro antialiasing analógico');
grid on;

%% FILTROS DIGITALES PARA 3 TASAS DE MUESTREO
filtros = {};
ordenes = [];
for i = 1:3
    fs = fs_list(i);
    fprintf('\n--- Filtro digital para fs = %d Hz ---\n', fs);
    Wp = fp / (fs/2);  % Normalización de la Frecuencia de paso de 1dB a la frecuencia de sampleo de cada señal
    [n, Wn] = buttord(Wp, 0.5, Ap, As);
    [b, a] = butter(n, Wn);
    filtros{i} = struct('fs', fs, 'b', b, 'a', a);
    ordenes(i) = n;
    % Graficar respuestas
    figure;
    freqz(b, a, 1024, fs);
    title(sprintf('Filtro digital Butterworth - fs = %d Hz', fs));
    drawnow; 

    figure;
    zplane(b, a);
    title(sprintf('Polos y ceros - fs = %d Hz', fs));
    drawnow; 
    
    figure;
    impz(b, a, [], fs);
    title(sprintf('Respuesta al impulso - fs = %d Hz', fs));
    drawnow; 
end

%% FILTRADO Y ESPECTRO DE SALIDA
fprintf('\n--- Filtrado digital de las señales ---\n')
for k = 1:3
    [x, fs] = audioread(archivos{k});
    fdata = filtros{k};
    y = filter(fdata.b, fdata.a, x); %Esta notacion toma los elementos del struct despues del punto

    % Espectro
    Y = abs(fft(y));
    N = length(Y);
    f = linspace(0, fs, N);

    figure(20);
    subplot(3,1,k);
    plot(f(1:N/2), 20*log10(Y(1:N/2)), 'Color', colores(k,:));
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(['Filtrada: ' archivos{k}]);
    grid on;
end


%% DECIMACIÓN Y CUANTIZACIÓN A 12 BITS DE TODAS LAS SEÑALES
fprintf('\n--- Decimación y cuantización para todas las señales ---\n')
for k = 1:3
    [x, fs] = audioread(archivos{k});
    fdata = filtros{k};
    y = filter(fdata.b, fdata.a, x);   % Aplicar el filtro digital

    % Decimación a 8 kHz
    factor_decimacion = fs / 8000;
    if mod(factor_decimacion,1) ~= 0
        error('La frecuencia de muestreo no es múltiplo de 8 kHz. No se puede decimar exactamente.');
    end
    y8k = decimate(y, factor_decimacion);

    % Cuantización a 12 bits
    y12b = round(y8k * (2^11 - 1));     % Rango [-2047, 2047]
    y12b = y12b / (2^11 - 1);           % Normalización [-1, 1]

    % Espectro
    Y = abs(fft(y12b));
    N = length(Y);
    f = linspace(0, 8000, N);

    figure(30);
    subplot(3,1,k);
    plot(f(1:N/2), 20*log10(Y(1:N/2)), 'Color', colores(k,:));
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(['Señal final (8kHz, 12 bits): ' archivos{k}]);
    grid on;

    % Reproducir
    fprintf('Reproduciendo original: %s\n', archivos{k});
    sound(x, fs); pause(3);
    fprintf('Reproduciendo procesada (8kHz): %s\n', archivos{k});
    sound(y12b, 8000); pause(3);
    % Guardar archivo procesado con sufijo '_procesada'
    [filepath, name, ext] = fileparts(archivos{k});
    nombre_salida = fullfile(filepath, [name '_procesada' ext]);
    audiowrite(nombre_salida, y12b, 8000);
    fprintf('Archivo guardado: %s\n', nombre_salida);
    
    % Archivos de entrada ya procesados a 8 kHz
archivos_procesados = {'16khz_procesada.wav', ...
                       '24khz_procesada.wav', ...
                       '32khz_procesada.wav'};

% Factor de interpolación
L = 15;
fs_in = 8000;
fs_out = fs_in * L;

    [s16p, fs] = audioread('16khz_procesada.wav'); 
    [s24p, fs] = audioread('24khz_procesada.wav'); 
    [s32p, fs] = audioread('32khz_procesada.wav'); 
    N  = length(s16p); %puede ser cualquiera porque miden lo mismo
    R = zeros(1, 15 * N);  % Vector de salida

for i = 1:N 

    posicion_muestra_no_nula = (i-1)*15 + 1;
    R(posicion_muestra_no_nula)   = s16p(i) + s24p(i) + s32p(i);

    
end
