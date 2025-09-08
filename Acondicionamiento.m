archivos16khz = {'canal1_16_khz.wav', 'canal2_16_khz.wav', 'canal3_16_khz.wav'}; %Incluyo arhivos nombrados por su fsampleo
archivos24khz = {'canal1_24_khz.wav', 'canal2_24_khz.wav', 'canal3_24_khz.wav'}; %Incluyo arhivos nombrados por su fsampleo
archivos32khz = {'canal1_32_khz.wav', 'canal2_32_khz.wav', 'canal3_32_khz.wav'}; %Incluyo arhivos nombrados por su fsampleo
archivos = {archivos16khz,archivos24khz,archivos32khz};
fs_list = [16000, 24000, 32000]; %Genero un vector con las fsamples (puedo hacerlo porque son conocidas)
colores = lines(3);  % Colores para graficar
load('filtros_guardados.mat', 'filtros');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
archivos_elegidos = archivos16khz;          
fdata = filtros{1}; % 16khz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DECIMACION Y CUANTIZACION A 12 BIT DE TODAS LAS SEÑALES
fprintf('\n--- Decimacion y cuantizacion para todas las señales ---\n')
for k = 1:3
    [x, fs] = audioread(archivos_elegidos{k});
    canalFiltrado = filter(fdata.b, fdata.a, x);   % Aplicar el filtro digital
    % Decimacion a 8 kHz
    factor_decimacion = fs / 8000;  %Se puede poner fijo elegido, pero por si cambia
    if mod(factor_decimacion,1) ~= 0
        error('La frecuencia de muestreo no es multiplo de 8 kHz. No se puede decimar exactamente.');
    end
    canalFiltrado8khz = decimate(canalFiltrado, factor_decimacion);

    % Cuantizacion a 12 bits
    canalFiltrado8khz12bits = round(canalFiltrado8khz * (2^11 - 1));     % Rango [-2047, 2047]
    canalFiltrado8khz12bits = canalFiltrado8khz12bits / (2^11 - 1);           % Normalizacion [-1, 1]

    % Espectro
    Y = abs(fft(canalFiltrado8khz12bits));
    N = length(Y);
    f = linspace(0,8000, N);

    figure(1);
    subplot(3,1,k);
    plot(f, 20*log10(Y), 'Color', colores(k,:));
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(['Canal Procesado (8kHz, 12 bits): ' archivos_elegidos{k}]);
    grid on;

    % Reproducir
    fprintf('Reproduciendo original: %s\n', archivos_elegidos{k});
    sound(x, fs); pause(3);
    fprintf('Reproduciendo procesada (8kHz): %s\n', archivos_elegidos{k});
    sound(canalFiltrado8khz12bits, 8000); pause(3);
    % Guardar archivo procesado con sufijo '_procesada'
    [filepath, name, ext] = fileparts(archivos_elegidos{k});
    nombre_salida = fullfile(filepath, ['canal_' int2str(k) ext]);
    audiowrite(nombre_salida, canalFiltrado8khz12bits, 8000);
    fprintf('Archivo guardado: %s\n', nombre_salida);    
    % Archivos de entrada ya procesados a 8 kHz
    archivos_procesados{k} = nombre_salida;
end

save('archivos_procesados.mat', 'archivos_procesados');