clear; clc;

%% CARGAR DATOS Y CONFIGURACIÓN
fprintf('=== DEMULTIPLEXOR FDM OPTIMIZADO - FILTRADO + DOWNSAMPLING EN UN PASO ===\n\n');

% Cargar la señal multiplexada y datos originales
load('salidaMux_polifasico.mat', 'salidaMux');
load('archivos_procesados.mat', 'archivos_procesados');
load('filtros_originales.mat', 'h_original');   % << Filtros originales

% Cargar señales originales para comparación
canal1_orig = audioread(archivos_procesados{1});
canal2_orig = audioread(archivos_procesados{2});
canal3_orig = audioread(archivos_procesados{3});

% Parámetros
fs_entrada_orig = 8000;   % Frecuencia original
fs_multiplexada = 120000; % Frecuencia de la señal multiplexada
factor = 15;              % Factor de downsampling

% Bandas de los canales
banda1 = [12300 15400];
banda2 = [16300 19400]; 
banda3 = [20300 23400];
bandas = {banda1, banda2, banda3};
nombres_canales = {'Canal 1', 'Canal 2', 'Canal 3'};
colores = {'b', 'r', 'g'};

L_mux = length(salidaMux);
L_orig = L_mux / factor;

fprintf('Longitud señal multiplexada: %d muestras\n', L_mux);
fprintf('Longitud esperada señales originales: %d muestras\n', L_orig);
fprintf('Factor de downsampling: %d\n\n', factor);

%% USAR FILTROS ORIGINALES
fprintf('--- Usando filtros originales guardados ---\n');
h_filtros = h_original;

for k = 1:3
    fprintf('Filtro %d cargado - Longitud: %d\n', k, length(h_filtros{k}));
end

%% ANÁLISIS ESPECTRAL DE LA SEÑAL MULTIPLEXADA
fprintf('\n--- Análisis de la señal multiplexada ---\n');

N_mux = length(salidaMux);
f_mux = linspace(0, fs_multiplexada, N_mux);
f_mux_plot = f_mux(1:N_mux/2);

X_mux = fft(salidaMux);
magX_mux_dB = 20*log10(abs(X_mux(1:N_mux/2)) + eps);

% Crear vector de tiempo para la señal multiplexada
t_mux = (0:N_mux-1)/fs_multiplexada;

%% GRÁFICA 1: SEÑAL MULTIPLEXADA DE ENTRADA
figure('Position', [100, 100, 1400, 600]);

subplot(1,2,1);
plot(t_mux(1:min(5000, end)), salidaMux(1:min(5000, end)), 'k', 'LineWidth', 1);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Multiplexada FDM - Entrada al Demultiplexor');
grid on;

subplot(1,2,2);
plot(f_mux_plot/1000, magX_mux_dB, 'k', 'LineWidth', 1.5);
xlabel('Frecuencia (kHz)');
ylabel('Magnitud (dB)');
title('Espectro de la Señal Multiplexada (120 kHz)');
grid on;
xlim([0, 60]);

% Resaltar las bandas de interés
hold on;
for k = 1:3
    ylims = ylim;
    fill([bandas{k}(1)/1000 bandas{k}(2)/1000 bandas{k}(2)/1000 bandas{k}(1)/1000], ...
         [ylims(1) ylims(1) ylims(2) ylims(2)], colores{k}, ...
         'FaceAlpha', 0.2, 'EdgeColor', colores{k}, 'LineStyle', '--');
    
    % Etiquetas
    text((bandas{k}(1)+bandas{k}(2))/2000, ylims(2)-5, ...
         sprintf('C%d', k), 'HorizontalAlignment', 'center', ...
         'FontWeight', 'bold', 'Color', colores{k});
end
hold off;

%% DEMULTIPLEXACIÓN OPTIMIZADA - FILTRADO + DOWNSAMPLING EN UN PASO
fprintf('\n--- Demultiplexación optimizada: Filtrado + Downsampling combinado ---\n');

% Inicializar señales demultiplexadas
canales_filtrados = cell(3,1);  % Para análisis posterior
canales_demux = cell(3,1);

% Procesar cada canal con método optimizado inline
for canal_idx = 1:3
    fprintf('Procesando %s (banda: %.1f-%.1f kHz)...\n', ...
            nombres_canales{canal_idx}, bandas{canal_idx}(1)/1000, bandas{canal_idx}(2)/1000);
    
    % MÉTODO OPTIMIZADO INLINE: Filtrado + Downsampling en un paso
    fprintf('  Aplicando filtrado + downsampling optimizado...\n');
    tic;
    
    % Obtener filtro actual
    h_actual = h_filtros{canal_idx};
    L_h = length(h_actual);
    L_x = length(salidaMux);
    L_out = floor(L_x / factor);
    
    % Inicializar señal de salida
    canal_recuperado = zeros(L_out, 1);
    
    fprintf('    Calculando %d muestras de salida (de %d posibles)...\n', L_out, L_x);
    
    % BUCLE PRINCIPAL OPTIMIZADO - Solo calcular muestras que se conservan
    for n = 1:L_out
        % Índice en la señal original correspondiente a esta muestra de salida
        n_original = (n-1) * factor + 1;
        
        % Calcular la convolución solo para este punto
        suma = 0;
        for k = 1:L_h
            idx_x = n_original - k + 1;
            if idx_x >= 1 && idx_x <= L_x
                suma = suma + h_actual(k) * salidaMux(idx_x);
            end
        end
        canal_recuperado(n) = suma;
        
        % Mostrar progreso cada 10000 muestras
        if mod(n, 10000) == 0 || n == L_out
            fprintf('      Progreso: %d/%d muestras (%.1f%%)\n', ...
                    n, L_out, (n/L_out)*100);
        end
    end
    
    tiempo_optimizado = toc;
    fprintf('  Tiempo de procesamiento: %.3f segundos\n', tiempo_optimizado);
    
    % APLICAR CORRECCIÓN DE DESPLAZAMIENTO EN FRECUENCIA (si es necesario)
    if canal_idx == 1 || canal_idx == 3
        fprintf('  Aplicando corrección de desplazamiento en frecuencia...\n');
        multDesplazamiento = 1;
        for n = 1:length(canal_recuperado)
            canal_recuperado(n) = canal_recuperado(n) * multDesplazamiento;
            multDesplazamiento = multDesplazamiento * (-1);
        end
    end
    
    canales_demux{canal_idx} = canal_recuperado;
    
    % Para análisis posterior, calcular señal filtrada completa (solo para gráficos)
    fprintf('  Calculando señal filtrada para análisis espectral...\n');
    senal_filtrada_completa = filter(h_actual, 1, salidaMux);
    canales_filtrados{canal_idx} = senal_filtrada_completa;
    
    fprintf('  Longitud canal recuperado: %d muestras\n', length(canal_recuperado));
    fprintf('  Longitud esperada: %d muestras\n\n', length(canal1_orig));
end

%% GUARDAR ARCHIVOS DE AUDIO RECUPERADOS
fprintf('--- Guardando archivos de audio recuperados ---\n');
for k = 1:3
    % Ajustar nombre de archivo
    nombre_archivo = sprintf('canal_%d_demux.wav', k);
    
    % Asegurar que los datos estén en rango [-1, 1]
    canal = canales_demux{k};
    
    % Normalización segura
    max_val = max(abs(canal));
    if max_val > 0
        canal = canal / max_val;
    end
    
    % Guardar audio
    audiowrite(nombre_archivo, canal, fs_entrada_orig);
    fprintf('Canal %d guardado en %s\n', k, nombre_archivo);
end

%% ANÁLISIS ESPECTRAL - PASO A PASO
fprintf('\n--- Análisis espectral paso a paso ---\n');

% Parámetros para FFT
N_orig = length(canales_demux{1});
f_orig = linspace(0, fs_entrada_orig, N_orig);
f_orig_plot = f_orig(1:N_orig/2);
t_orig = (0:N_orig-1)/fs_entrada_orig;

%% GRÁFICA 2: PROCESO PASO A PASO PARA CADA CANAL
figure('Position', [200, 200, 1600, 1200]);

canales_originales = {canal1_orig, canal2_orig, canal3_orig};

for k = 1:3
    % PASO 1: Respuesta del filtro
    subplot(3,4,4*k-3);
    [H_resp, w_resp] = freqz(h_filtros{k}, 1, 2048, fs_multiplexada);
    plot(w_resp/1000, 20*log10(abs(H_resp)), colores{k}, 'LineWidth', 2);
    xlabel('Frecuencia (kHz)');
    ylabel('Magnitud (dB)');
    title(sprintf('%s - Respuesta Filtro', nombres_canales{k}));
    grid on;
    xlim([0, 30]);
    ylim([-80, 10]);
    
    % Resaltar banda objetivo
    hold on;
    ylims = ylim;
    fill([bandas{k}(1)/1000 bandas{k}(2)/1000 bandas{k}(2)/1000 bandas{k}(1)/1000], ...
         [ylims(1) ylims(1) ylims(2) ylims(2)], colores{k}, ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold off;
    
    % PASO 2: Señal filtrada (120 kHz)
    subplot(3,4,4*k-2);
    X_filt = fft(canales_filtrados{k});
    magX_filt_dB = 20*log10(abs(X_filt(1:length(X_filt)/2)) + eps);
    f_filt = linspace(0, fs_multiplexada, length(X_filt));
    f_filt_plot = f_filt(1:length(f_filt)/2);
    
    plot(f_filt_plot/1000, magX_filt_dB, colores{k}, 'LineWidth', 1.5);
    xlabel('Frecuencia (kHz)');
    ylabel('Magnitud (dB)');
    title('Después Filtrado (120 kHz)');
    grid on;
    xlim([0, 30]);
    
    % PASO 3: Después del downsampling optimizado
    subplot(3,4,4*k-1);
    X_demux = fft(canales_demux{k});
    magX_demux_dB = 20*log10(abs(X_demux(1:N_orig/2)) + eps);
    
    plot(f_orig_plot, magX_demux_dB, colores{k}, 'LineWidth', 2);
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title('Método Optimizado (8 kHz)');
    grid on;
    xlim([0, 4000]);
    
    % PASO 4: Comparación con original
    subplot(3,4,4*k);
    X_orig = fft(canales_originales{k}(1:N_orig));
    magX_orig_dB = 20*log10(abs(X_orig(1:N_orig/2)) + eps);
    
    plot(f_orig_plot, magX_orig_dB, 'k', 'LineWidth', 2, 'DisplayName', 'Original');
    hold on;
    plot(f_orig_plot, magX_demux_dB, colores{k}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Recuperado');
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title('Comparación Final');
    legend('show');
    grid on;
    xlim([0, 4000]);
    hold off;
end

%% GRÁFICA 3: COMPARACIÓN TEMPORAL
figure('Position', [300, 300, 1400, 900]);

for k = 1:3
    % Ajustar longitudes para comparación
    len_min = min(length(canales_demux{k}), length(canales_originales{k}));
    t_comp = (0:len_min-1)/fs_entrada_orig;
    t_plot = t_comp(1:min(5000, len_min));
    
    subplot(3,2,2*k-1);
    plot(t_plot, canales_originales{k}(1:length(t_plot)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Original');
    hold on;
    plot(t_plot, canales_demux{k}(1:length(t_plot)), colores{k}, 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Recuperado');
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    title([nombres_canales{k} ' - Comparación Temporal']);
    legend('show');
    grid on;
    hold off;
    
    % Error temporal
    subplot(3,2,2*k);
    error_temporal = canales_demux{k}(1:len_min) - canales_originales{k}(1:len_min);
    plot(t_comp(1:min(5000, len_min)), error_temporal(1:min(5000, len_min)), 'k', 'LineWidth', 1);
    xlabel('Tiempo (s)');
    ylabel('Error');
    title([nombres_canales{k} ' - Error de Recuperación']);
    grid on;
    
    % Estadísticas en el gráfico
    error_rms = sqrt(mean(error_temporal.^2));
    error_max = max(abs(error_temporal));
    text(0.02, 0.95, sprintf('RMS: %.4f\nMáx: %.4f', error_rms, error_max), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', 'BackgroundColor', 'white');
end

%% ANÁLISIS CUANTITATIVO DETALLADO
fprintf('\n=== ANÁLISIS CUANTITATIVO DE RECUPERACIÓN ===\n');

correlaciones = zeros(3,1);
snr_valores = zeros(3,1);
errores_relativos = zeros(3,1);

for k = 1:3
    fprintf('\n%s:\n', nombres_canales{k});
    fprintf('----------------------------------------\n');
    
    % Ajustar longitudes
    len_min = min(length(canales_demux{k}), length(canales_originales{k}));
    canal_demux_ajustado = canales_demux{k}(1:len_min);
    canal_orig_ajustado = canales_originales{k}(1:len_min);
    
    % Errores temporales
    error_temporal = canal_demux_ajustado - canal_orig_ajustado;
    error_rms = sqrt(mean(error_temporal.^2));
    error_max = max(abs(error_temporal));
    
    % Potencias
    potencia_orig = sqrt(mean(canal_orig_ajustado.^2));
    potencia_demux = sqrt(mean(canal_demux_ajustado.^2));
    
    % Error relativo y SNR
    error_relativo = error_rms / potencia_orig;
    snr_db = 20*log10(potencia_orig / error_rms);
    
    % Correlación
    correlacion = corrcoef(canal_demux_ajustado, canal_orig_ajustado);
    correlacion_valor = correlacion(1,2);
    
    % Guardar para estadísticas globales
    correlaciones(k) = correlacion_valor;
    snr_valores(k) = snr_db;
    errores_relativos(k) = error_relativo;
    
    fprintf('  ERRORES TEMPORALES:\n');
    fprintf('    Error RMS: %.6f\n', error_rms);
    fprintf('    Error máximo: %.6f\n', error_max);
    fprintf('    Error relativo: %.4f%%\n', error_relativo * 100);
    
    fprintf('  MÉTRICAS DE CALIDAD:\n');
    fprintf('    Correlación: %.6f\n', correlacion_valor);
    fprintf('    SNR: %.2f dB\n', snr_db);
    fprintf('    Potencia original: %.6f\n', potencia_orig);
    fprintf('    Potencia recuperada: %.6f\n', potencia_demux);
    fprintf('    Relación potencias: %.4f\n', potencia_demux/potencia_orig);
    
    % Análisis espectral
    X_orig_full = fft(canal_orig_ajustado);
    X_demux_full = fft(canal_demux_ajustado);
    
    % Energía total
    energia_orig = sum(abs(X_orig_full).^2);
    energia_demux = sum(abs(X_demux_full).^2);
    
    % Energía en banda de audio (300-3400 Hz)
    f_analisis = linspace(0, fs_entrada_orig, len_min);
    idx_audio = find(f_analisis >= 300 & f_analisis <= 3400);
    energia_orig_audio = sum(abs(X_orig_full(idx_audio)).^2);
    energia_demux_audio = sum(abs(X_demux_full(idx_audio)).^2);
    
    eficiencia_total = energia_demux / energia_orig;
    eficiencia_audio = energia_demux_audio / energia_orig_audio;
    
    fprintf('  ANÁLISIS ESPECTRAL:\n');
    fprintf('    Eficiencia total: %.4f (%.2f%%)\n', eficiencia_total, eficiencia_total*100);
    fprintf('    Eficiencia banda audio: %.4f (%.2f%%)\n', eficiencia_audio, eficiencia_audio*100);
    
    % Evaluar calidad general
    if correlacion_valor > 0.95 && error_relativo < 0.05
        calidad = 'EXCELENTE';
    elseif correlacion_valor > 0.90 && error_relativo < 0.10
        calidad = 'BUENA';
    elseif correlacion_valor > 0.80 && error_relativo < 0.20
        calidad = 'ACEPTABLE';
    else
        calidad = 'DEFICIENTE';
    end
    
    fprintf('  EVALUACIÓN GENERAL: %s\n', calidad);
end

%% GRÁFICA 4: ANÁLISIS DE ALIASING Y EFECTOS
figure('Position', [400, 400, 1400, 800]);

for k = 1:3
    subplot(2,3,k);
    
    % Espectro del canal filtrado a 120 kHz (antes del downsampling)
    X_filt = fft(canales_filtrados{k});
    magX_filt_dB = 20*log10(abs(X_filt(1:length(X_filt)/2)) + eps);
    f_filt = linspace(0, fs_multiplexada, length(X_filt));
    f_filt_plot = f_filt(1:length(f_filt)/2);
    
    plot(f_filt_plot/1000, magX_filt_dB, colores{k}, 'LineWidth', 1.5);
    xlabel('Frecuencia (kHz)');
    ylabel('Magnitud (dB)');
    title([nombres_canales{k} ' - Antes Downsampling']);
    grid on;
    xlim([0, 60]);
    
    % Marcar frecuencias de aliasing potencial
    hold on;
    for aliasing_freq = fs_entrada_orig/1000:fs_entrada_orig/1000:60
        ylims_current = ylim;
        line([aliasing_freq aliasing_freq], ylims_current, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 0.5);
    end
    
    % Marcar banda original
    band_center = (bandas{k}(1) + bandas{k}(2))/2000;
    ylims_current = ylim;
    line([band_center band_center], ylims_current, 'Color', colores{k}, 'LineStyle', '--', 'LineWidth', 2);
    hold off;
    
    % Análisis después del downsampling
    subplot(2,3,k+3);
    len_min = min(length(canales_demux{k}), length(canales_originales{k}));
    X_demux = fft(canales_demux{k}(1:len_min));
    X_orig = fft(canales_originales{k}(1:len_min));
    
    magX_demux_dB = 20*log10(abs(X_demux(1:len_min/2)) + eps);
    magX_orig_dB = 20*log10(abs(X_orig(1:len_min/2)) + eps);
    f_comp = linspace(0, fs_entrada_orig, len_min);
    f_comp_plot = f_comp(1:len_min/2);
    
    plot(f_comp_plot, magX_orig_dB, colores{k}, 'LineWidth', 2, 'DisplayName', 'Original');
    hold on;
    plot(f_comp_plot, magX_demux_dB, colores{k}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Recuperado');
    
    % Mostrar diferencia
    diferencia_dB = magX_demux_dB - magX_orig_dB;
    plot(f_comp_plot, diferencia_dB - 20, 'r:', 'LineWidth', 1, 'DisplayName', 'Diferencia-20dB');
    
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title([nombres_canales{k} ' - Comparación Final']);
    legend('show');
    grid on;
    xlim([0, 4000]);
    hold off;
end

%% RESUMEN FINAL
fprintf('\n=== RESUMEN DE LA DEMULTIPLEXACIÓN OPTIMIZADA ===\n');
fprintf('✓ Método: Filtrado + Downsampling combinado en un solo paso\n');
fprintf('✓ Ventajas del método optimizado:\n');
fprintf('  - Menor uso de memoria\n');
fprintf('  - Mayor eficiencia computacional\n');
fprintf('  - Solo calcula las muestras que se van a conservar\n');
fprintf('✓ %d canales procesados exitosamente\n', length(canales_demux));

fprintf('\nESTADÍSTICAS GLOBALES:\n');
fprintf('  Correlación promedio: %.4f\n', mean(correlaciones));
fprintf('  SNR promedio: %.2f dB\n', mean(snr_valores));
fprintf('  Error relativo promedio: %.4f%%\n', mean(errores_relativos)*100);

if mean(correlaciones) > 0.9
    fprintf('  ✅ RECUPERACIÓN EXITOSA - Las señales se recuperan satisfactoriamente\n');
elseif mean(correlaciones) > 0.7
    fprintf('  ⚠️ RECUPERACIÓN PARCIAL - Las señales se recuperan con cierta degradación\n');
else
    fprintf('  ❌ RECUPERACIÓN DEFICIENTE - Las señales presentan alta distorsión\n');
end

%% GUARDAR RESULTADOS
save('resultados_demultiplexacion_optimizada.mat', 'canales_demux', 'canales_filtrados', 'correlaciones', 'snr_valores', 'errores_relativos');

fprintf('\n✓ Resultados guardados en: resultados_demultiplexacion_optimizada.mat\n');
fprintf('✓ Archivos de audio guardados como: canal_1_demux.wav, canal_2_demux.wav, canal_3_demux.wav\n');
fprintf('✓ Análisis completado con método optimizado\n\n');