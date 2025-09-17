


% Par�metros
fs_salida = 120000; % Frecuencia de muestreo de la se�al multiplexada

% Definir las bandas de inter�s para la visualizaci�n
banda1 = [12300 15400];
banda2 = [16300 19400];
banda3 = [20300 23400];
bandas = {banda1, banda2, banda3};
nombres_canales = {'Canal 1 (12.3-15.4 kHz)', 'Canal 2 (16.3-19.4 kHz)', 'Canal 3 (20.3-23.4 kHz)'};
colores = {'b', 'r', 'g'};

%% Aplicar filtros a la se�al multiplexada
% Inicializar un cell array para guardar las salidas
canales_demultiplexados = cell(3,1);

for k = 1:3
    fprintf('Aplicando filtro del %s a la se�al multiplexada...\n', nombres_canales{k});
    
    % Aplicar el filtro original a la se�al multiplexada
    canales_demultiplexados{k} = filter(h_original{k}, 1, salidaMux);
end

%% An�lisis y visualizaci�n de espectros
N = length(salidaMux);
f = linspace(0, fs_salida, N);
f_nyquist = f(1:N/2);

figure('Position', [100, 100, 1200, 800]);

for k = 1:3
    % Calcular el espectro de la se�al demultiplexada
    X_demux = fft(canales_demultiplexados{k});
    magX_demux_dB = 20*log10(abs(X_demux(1:N/2)) + eps);
    
    % Subgr�fica para cada canal
    subplot(3,1,k);
    plot(f_nyquist/1000, magX_demux_dB, colores{k}, 'LineWidth', 1.5);
    
    % Resaltar la banda de inter�s
    hold on;
    ylims = ylim;
    fill([bandas{k}(1)/1000, bandas{k}(2)/1000, bandas{k}(2)/1000, bandas{k}(1)/1000], ...
         [ylims(1), ylims(1), ylims(2), ylims(2)], colores{k}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
    xlabel('Frecuencia (kHz)');
    ylabel('Magnitud (dB)');
    title(sprintf('Espectro de la se�al recuperada para el %s', nombres_canales{k}));
    grid on;
    xlim([0, 30]); % Limitar la visualizaci�n hasta 30 kHz
    hold off;
end



%% Gr�fica superpuesta para comparar
figure;
hold on;
for k = 1:3
    X_demux = fft(canales_demultiplexados{k});
    magX_demux_dB = 20*log10(abs(X_demux(1:N/2)) + eps);
    plot(f_nyquist/1000, magX_demux_dB, colores{k}, 'LineWidth', 1.5, 'DisplayName', nombres_canales{k});
end

xlabel('Frecuencia (kHz)');
ylabel('Magnitud (dB)');
title('Comparaci�n de espectros de canales demultiplexados');
legend('show', 'Location', 'best');
grid on;
xlim([0, 30]); % Limitar la visualizaci�n
hold off;

fprintf('\nProceso de demultiplexado y an�lisis espectral completado.\n');









% Cargar las se�ales demultiplexadas si no est�n en el workspace
if ~exist('canales_demultiplexados', 'var')
    load('canales_demultiplexados.mat', 'canales_demultiplexados');
    if ~exist('canales_demultiplexados', 'var')
        error('Variable canales_demultiplexados no encontrada. Ejecute el script anterior primero.');
    end
end

% Par�metros de downsampling
factor_decimacion = 15;
fs_salida = 120000;
fs_downsampled = fs_salida / factor_decimacion; % 120000 / 15 = 8000 Hz

nombres_canales = {'Canal 1 (12.3-15.4 kHz)', 'Canal 2 (16.3-19.4 kHz)', 'Canal 3 (20.3-23.4 kHz)'};
colores = {'b', 'r', 'g'};

%% Downsampling de cada se�al
canales_downsampled = cell(3, 1);

for k = 1:3
    fprintf('Realizando downsampling del %s (factor %d)...\n', nombres_canales{k}, factor_decimacion);
    
    % Downsampling simple tomando 1 de cada 15 muestras
    % OJO: Esta es una decimaci�n simple sin filtro anti-aliasing.
    canales_downsampled{k} = canales_demultiplexados{k}(1:factor_decimacion:end);
end

%% Desplazar la frecuencia de los canales 1 y 3 de vuelta a banda base
% El downsampling no cambia la posici�n de la frecuencia.
% El espectro se "pliega" (aliasing) sobre s� mismo. Para recuperar la
% se�al original, se debe aplicar el mismo desplazamiento de frecuencia
% inverso que se aplic� en el multiplexado.
multDesplazamintoFrecuencia = 1;
for k = [1, 3] % Solo canales 1 y 3
    for n = 1:length(canales_downsampled{k})
        canales_downsampled{k}(n) = canales_downsampled{k}(n) * multDesplazamintoFrecuencia;
        multDesplazamintoFrecuencia = multDesplazamintoFrecuencia * (-1);
    end
    multDesplazamintoFrecuencia = 1; % Reiniciar para el siguiente canal
end

%% An�lisis y visualizaci�n de espectros
N = length(canales_downsampled{1});
f = linspace(0, fs_downsampled, N);
f_nyquist_ds = f(1:N/2);

figure('Position', [100, 100, 1200, 800]);

for k = 1:3
    % Calcular el espectro de la se�al downsampleada
    X_ds = fft(canales_downsampled{k});
    magX_ds_dB = 20*log10(abs(X_ds(1:N/2)) + eps);
    
    % Subgr�fica para cada canal
    subplot(3,1,k);
    plot(f_nyquist_ds, magX_ds_dB, colores{k}, 'LineWidth', 1.5);
    
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(sprintf('Espectro del canal recuperado %d (fs = 8 kHz)', k));
    grid on;
    xlim([0, fs_downsampled/2]); % Mostrar hasta la nueva frecuencia de Nyquist (4 kHz)
end



%% Comparaci�n superpuesta
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
title('Comparaci�n de espectros de los 3 canales recuperados');
legend('show', 'Location', 'best');
grid on;
xlim([0, fs_downsampled/2]);

fprintf('\nProceso de downsampling y an�lisis espectral de las se�ales recuperadas completado.\n');

% Cargar las se�ales si no est�n en el workspace
if ~exist('canales_downsampled', 'var')
    load('canales_recuperados.mat', 'canales_downsampled'); 
    if ~exist('canales_downsampled', 'var')
        error('Variable canales_downsampled no encontrada. Aseg�rese de haber ejecutado los scripts anteriores.');
    end
end

% Par�metros de guardado
fs_audio = 8000; % Frecuencia de muestreo de las se�ales recuperadas
nombres_archivos = {'canal1_recuperado.wav', 'canal2_recuperado.wav', 'canal3_recuperado.wav'};
nombres_canales = {'Canal 1', 'Canal 2', 'Canal 3'};

fprintf('--- Guardando archivos de audio --- \n');

for k = 1:3
    % Obtener la se�al del canal
    canal_recuperado = canales_downsampled{k};
    
    % Normalizar la se�al para evitar distorsi�n
    % La normalizaci�n escala los valores entre -1 y 1
    canal_normalizado = canal_recuperado / max(abs(canal_recuperado));
    
    % Guardar la se�al en un archivo .wav
    try
        audiowrite(nombres_archivos{k}, canal_normalizado, fs_audio);
        fprintf('Se�al del %s guardada exitosamente como: %s\n', nombres_canales{k}, nombres_archivos{k});
    catch ME
        fprintf('Error al guardar el archivo %s: %s\n', nombres_archivos{k}, ME.message);
    end
end

fprintf('\nProceso completado. �Ahora puedes reproducir los archivos de audio!\n');