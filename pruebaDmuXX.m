


% Parámetros
fs_salida = 120000; % Frecuencia de muestreo de la señal multiplexada

% Definir las bandas de interés para la visualización
banda1 = [12300 15400];
banda2 = [16300 19400];
banda3 = [20300 23400];
bandas = {banda1, banda2, banda3};
nombres_canales = {'Canal 1 (12.3-15.4 kHz)', 'Canal 2 (16.3-19.4 kHz)', 'Canal 3 (20.3-23.4 kHz)'};
colores = {'b', 'r', 'g'};

%% Aplicar filtros a la señal multiplexada
% Inicializar un cell array para guardar las salidas
canales_demultiplexados = cell(3,1);

for k = 1:3
    fprintf('Aplicando filtro del %s a la señal multiplexada...\n', nombres_canales{k});
    
    % Aplicar el filtro original a la señal multiplexada
    canales_demultiplexados{k} = filter(h_original{k}, 1, salidaMux);
end

%% Análisis y visualización de espectros
N = length(salidaMux);
f = linspace(0, fs_salida, N);
f_nyquist = f(1:N/2);

figure('Position', [100, 100, 1200, 800]);

for k = 1:3
    % Calcular el espectro de la señal demultiplexada
    X_demux = fft(canales_demultiplexados{k});
    magX_demux_dB = 20*log10(abs(X_demux(1:N/2)) + eps);
    
    % Subgráfica para cada canal
    subplot(3,1,k);
    plot(f_nyquist/1000, magX_demux_dB, colores{k}, 'LineWidth', 1.5);
    
    % Resaltar la banda de interés
    hold on;
    ylims = ylim;
    fill([bandas{k}(1)/1000, bandas{k}(2)/1000, bandas{k}(2)/1000, bandas{k}(1)/1000], ...
         [ylims(1), ylims(1), ylims(2), ylims(2)], colores{k}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
    xlabel('Frecuencia (kHz)');
    ylabel('Magnitud (dB)');
    title(sprintf('Espectro de la señal recuperada para el %s', nombres_canales{k}));
    grid on;
    xlim([0, 30]); % Limitar la visualización hasta 30 kHz
    hold off;
end



%% Gráfica superpuesta para comparar
figure;
hold on;
for k = 1:3
    X_demux = fft(canales_demultiplexados{k});
    magX_demux_dB = 20*log10(abs(X_demux(1:N/2)) + eps);
    plot(f_nyquist/1000, magX_demux_dB, colores{k}, 'LineWidth', 1.5, 'DisplayName', nombres_canales{k});
end

xlabel('Frecuencia (kHz)');
ylabel('Magnitud (dB)');
title('Comparación de espectros de canales demultiplexados');
legend('show', 'Location', 'best');
grid on;
xlim([0, 30]); % Limitar la visualización
hold off;

fprintf('\nProceso de demultiplexado y análisis espectral completado.\n');









% Cargar las señales demultiplexadas si no están en el workspace
if ~exist('canales_demultiplexados', 'var')
    load('canales_demultiplexados.mat', 'canales_demultiplexados');
    if ~exist('canales_demultiplexados', 'var')
        error('Variable canales_demultiplexados no encontrada. Ejecute el script anterior primero.');
    end
end

% Parámetros de downsampling
factor_decimacion = 15;
fs_salida = 120000;
fs_downsampled = fs_salida / factor_decimacion; % 120000 / 15 = 8000 Hz

nombres_canales = {'Canal 1 (12.3-15.4 kHz)', 'Canal 2 (16.3-19.4 kHz)', 'Canal 3 (20.3-23.4 kHz)'};
colores = {'b', 'r', 'g'};

%% Downsampling de cada señal
canales_downsampled = cell(3, 1);

for k = 1:3
    fprintf('Realizando downsampling del %s (factor %d)...\n', nombres_canales{k}, factor_decimacion);
    
    % Downsampling simple tomando 1 de cada 15 muestras
    % OJO: Esta es una decimación simple sin filtro anti-aliasing.
    canales_downsampled{k} = canales_demultiplexados{k}(1:factor_decimacion:end);
end

%% Desplazar la frecuencia de los canales 1 y 3 de vuelta a banda base
% El downsampling no cambia la posición de la frecuencia.
% El espectro se "pliega" (aliasing) sobre sí mismo. Para recuperar la
% señal original, se debe aplicar el mismo desplazamiento de frecuencia
% inverso que se aplicó en el multiplexado.
multDesplazamintoFrecuencia = 1;
for k = [1, 3] % Solo canales 1 y 3
    for n = 1:length(canales_downsampled{k})
        canales_downsampled{k}(n) = canales_downsampled{k}(n) * multDesplazamintoFrecuencia;
        multDesplazamintoFrecuencia = multDesplazamintoFrecuencia * (-1);
    end
    multDesplazamintoFrecuencia = 1; % Reiniciar para el siguiente canal
end

%% Análisis y visualización de espectros
N = length(canales_downsampled{1});
f = linspace(0, fs_downsampled, N);
f_nyquist_ds = f(1:N/2);

figure('Position', [100, 100, 1200, 800]);

for k = 1:3
    % Calcular el espectro de la señal downsampleada
    X_ds = fft(canales_downsampled{k});
    magX_ds_dB = 20*log10(abs(X_ds(1:N/2)) + eps);
    
    % Subgráfica para cada canal
    subplot(3,1,k);
    plot(f_nyquist_ds, magX_ds_dB, colores{k}, 'LineWidth', 1.5);
    
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(sprintf('Espectro del canal recuperado %d (fs = 8 kHz)', k));
    grid on;
    xlim([0, fs_downsampled/2]); % Mostrar hasta la nueva frecuencia de Nyquist (4 kHz)
end



%% Comparación superpuesta
figure;
hold on;
for k = 1:3
    X_ds = fft(canales_downsampled{k});
    magX_ds_dB = 20*log10(abs(X_ds(1:N/2)) + eps);
    plot(f_nyquist_ds, magX_ds_dB, colores{k}, 'LineWidth', 1.5, 'DisplayName', nombres_canales{k});
end
hold off;

xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Comparación de espectros de los 3 canales recuperados');
legend('show', 'Location', 'best');
grid on;
xlim([0, fs_downsampled/2]);

fprintf('\nProceso de downsampling y análisis espectral de las señales recuperadas completado.\n');

% Cargar las señales si no están en el workspace
if ~exist('canales_downsampled', 'var')
    load('canales_recuperados.mat', 'canales_downsampled'); 
    if ~exist('canales_downsampled', 'var')
        error('Variable canales_downsampled no encontrada. Asegúrese de haber ejecutado los scripts anteriores.');
    end
end

% Parámetros de guardado
fs_audio = 8000; % Frecuencia de muestreo de las señales recuperadas
nombres_archivos = {'canal1_recuperado.wav', 'canal2_recuperado.wav', 'canal3_recuperado.wav'};
nombres_canales = {'Canal 1', 'Canal 2', 'Canal 3'};

fprintf('--- Guardando archivos de audio --- \n');

for k = 1:3
    % Obtener la señal del canal
    canal_recuperado = canales_downsampled{k};
    
    % Normalizar la señal para evitar distorsión
    % La normalización escala los valores entre -1 y 1
    canal_normalizado = canal_recuperado / max(abs(canal_recuperado));
    
    % Guardar la señal en un archivo .wav
    try
        audiowrite(nombres_archivos{k}, canal_normalizado, fs_audio);
        fprintf('Señal del %s guardada exitosamente como: %s\n', nombres_canales{k}, nombres_archivos{k});
    catch ME
        fprintf('Error al guardar el archivo %s: %s\n', nombres_archivos{k}, ME.message);
    end
end

fprintf('\nProceso completado. ¡Ahora puedes reproducir los archivos de audio!\n');