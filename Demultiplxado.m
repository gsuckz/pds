clear; clc;

%% CARGAR DATOS Y CONFIGURACIÓN
fprintf('=== DEMULTIPLEXOR FDM OPTIMIZADO - FILTRADO + DOWNSAMPLING EN UN PASO ===\n\n');

% Cargar la señal multiplexada y datos originales
load('salidaMux_polifasico.mat', 'salidaMux');
load('archivos_procesados.mat', 'archivos_procesados');
load('filtros_originales.mat', 'h_original');

% Cargar señales originales para comparación (UNA SOLA VEZ)
canales_originales = cell(3,1);
canales_originales{1} = audioread(archivos_procesados{1});
canales_originales{2} = audioread(archivos_procesados{2});
canales_originales{3} = audioread(archivos_procesados{3});

% Parámetros
fs_entrada_orig = 8000;
fs_multiplexada = 120000;
factor = 15;

% Bandas de los canales
bandas = {[12300 15400], [16300 19400], [20300 23400]};
nombres_canales = {'Canal 1', 'Canal 2', 'Canal 3'};
colores = {'b', 'r', 'g'};

% Precálculo de constantes
L_mux = length(salidaMux);
L_orig = L_mux / factor;
N_mux = L_mux;
N_orig = floor(L_orig);

fprintf('Longitud señal multiplexada: %d muestras\n', L_mux);
fprintf('Longitud esperada señales originales: %d muestras\n', N_orig);
fprintf('Factor de downsampling: %d\n\n', factor);

%% USAR FILTROS ORIGINALES
fprintf('--- Usando filtros originales guardados ---\n');
h_filtros = h_original;
for k = 1:3
    fprintf('Filtro %d cargado - Longitud: %d\n', k, length(h_filtros{k}));
end

%% ANÁLISIS ESPECTRAL DE LA SEÑAL MULTIPLEXADA
fprintf('\n--- Análisis de la señal multiplexada ---\n');

% Precalcular vectores de frecuencia (UNA SOLA VEZ)
f_mux = linspace(0, fs_multiplexada, N_mux);
f_mux_plot = f_mux(1:N_mux/2);

% FFT de entrada (calcular UNA SOLA VEZ)
X_mux = fft(salidaMux);
magX_mux_dB = 20*log10(abs(X_mux(1:N_mux/2)) + eps);

t_mux = (0:N_mux-1)/fs_multiplexada;

%% GRÁFICA 1: SEÑAL MULTIPLEXADA DE ENTRADA
figure('Position', [100, 100, 1400, 600]);

subplot(1,2,1);
idx_plot = 1:min(5000, N_mux);
plot(t_mux(idx_plot), salidaMux(idx_plot), 'k', 'LineWidth', 1);
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

hold on;
ylims = ylim;
for k = 1:3
    fill([bandas{k}(1)/1000 bandas{k}(2)/1000 bandas{k}(2)/1000 bandas{k}(1)/1000], ...
         [ylims(1) ylims(1) ylims(2) ylims(2)], colores{k}, ...
         'FaceAlpha', 0.2, 'EdgeColor', colores{k}, 'LineStyle', '--');
    text((bandas{k}(1)+bandas{k}(2))/2000, ylims(2)-5, ...
         sprintf('C%d', k), 'HorizontalAlignment', 'center', ...
         'FontWeight', 'bold', 'Color', colores{k});
end
hold off;

%% DEMULTIPLEXACIÓN OPTIMIZADA
fprintf('\n--- Demultiplexación optimizada: Filtrado + Downsampling combinado ---\n');

canales_demux = cell(3,1);
canales_filtrados = cell(3,1);

% Precalcular L_out (es igual para todos los canales)
L_out = floor(L_mux / factor);

for canal_idx = 1:3
    fprintf('Procesando %s (banda: %.1f-%.1f kHz)...\n', ...
            nombres_canales{canal_idx}, bandas{canal_idx}(1)/1000, bandas{canal_idx}(2)/1000);
    
    fprintf('  Aplicando filtrado + downsampling optimizado...\n');
    tic;
    
    % Obtener filtro actual
    h_actual = h_filtros{canal_idx};
    L_h = length(h_actual);
    
    % OPTIMIZACIÓN CRÍTICA: Pre-alocar memoria
    canal_recuperado = zeros(L_out, 1);
    
    % BUCLE OPTIMIZADO - Menos accesos a índices
    fprintf('    Calculando %d muestras de salida...\n', L_out);
    
    % Optimización: Calcular convulación directa eficiente
    for n = 1:L_out
        n_original = (n-1) * factor + 1;
        
        suma = 0;
        for k = 1:L_h
            idx_x = n_original - k + 1;
            if idx_x >= 1 && idx_x <= L_mux
                suma = suma + h_actual(k) * salidaMux(idx_x);
            end
        end
        canal_recuperado(n) = suma;
        
        % Progreso reducido para menos I/O
        if mod(n, 20000) == 0 || n == L_out
            fprintf('      Progreso: %d/%d (%.1f%%)\n', n, L_out, (n/L_out)*100);
        end
    end
    
    tiempo_opt = toc;
    fprintf('  Tiempo: %.3f segundos\n', tiempo_opt);
    
    % CORRECCIÓN DE DESPLAZAMIENTO optimizada
    if canal_idx == 1 || canal_idx == 3
        fprintf('  Aplicando corrección de desplazamiento...\n');
        % Multiplicación vectorizada (MÁS EFICIENTE)
        canal_recuperado = canal_recuperado .* ((-1).^(0:L_out-1)');
    end
    
    canales_demux{canal_idx} = canal_recuperado;
    
    % Calcular señal filtrada SOLO para análisis (usando filter optimizado de MATLAB)
    fprintf('  Calculando señal filtrada para análisis...\n');
    canales_filtrados{canal_idx} = filter(h_actual, 1, salidaMux);
    
    fprintf('  Longitud recuperada: %d muestras\n\n', length(canal_recuperado));
end

%% GUARDAR ARCHIVOS DE AUDIO
fprintf('--- Guardando archivos de audio recuperados ---\n');
for k = 1:3
    nombre_archivo = sprintf('canal_%d_demux.wav', k);
    canal = canales_demux{k};
    
    % Normalización vectorizada
    max_val = max(abs(canal));
    if max_val > 0
        canal = canal / max_val;
    end
    
    audiowrite(nombre_archivo, canal, fs_entrada_orig);
    fprintf('Canal %d guardado en %s\n', k, nombre_archivo);
end

%% ANÁLISIS ESPECTRAL - OPTIMIZADO
fprintf('\n--- Análisis espectral paso a paso ---\n');

% Precalcular vectores de frecuencia para 8 kHz (UNA SOLA VEZ)
f_orig = linspace(0, fs_entrada_orig, N_orig);
f_orig_plot = f_orig(1:N_orig/2);
t_orig = (0:N_orig-1)/fs_entrada_orig;

%% GRÁFICA 2: PROCESO PASO A PASO
figure('Position', [200, 200, 1600, 1200]);

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
    
    hold on;
    ylims = ylim;
    fill([bandas{k}(1)/1000 bandas{k}(2)/1000 bandas{k}(2)/1000 bandas{k}(1)/1000], ...
         [ylims(1) ylims(1) ylims(2) ylims(2)], colores{k}, ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold off;
    
    % PASO 2: Señal filtrada
    subplot(3,4,4*k-2);
    X_filt = fft(canales_filtrados{k});
    L_filt = length(X_filt);
    magX_filt_dB = 20*log10(abs(X_filt(1:L_filt/2)) + eps);
    f_filt_plot = linspace(0, fs_multiplexada/2, L_filt/2);
    
    plot(f_filt_plot/1000, magX_filt_dB, colores{k}, 'LineWidth', 1.5);
    xlabel('Frecuencia (kHz)');
    ylabel('Magnitud (dB)');
    title('Después Filtrado (120 kHz)');
    grid on;
    xlim([0, 30]);
    
    % PASO 3: Después del downsampling
    subplot(3,4,4*k-1);
    X_demux = fft(canales_demux{k});
    magX_demux_dB = 20*log10(abs(X_demux(1:N_orig/2)) + eps);
    
    plot(f_orig_plot, magX_demux_dB, colores{k}, 'LineWidth', 2);
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title('Método Optimizado (8 kHz)');
    grid on;
    xlim([0, 4000]);
    
    % PASO 4: Comparación
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

% Precalcular longitud mínima común
len_min = min([length(canales_demux{1}), length(canales_demux{2}), length(canales_demux{3}), ...
               length(canales_originales{1}), length(canales_originales{2}), length(canales_originales{3})]);
t_comp = (0:len_min-1)/fs_entrada_orig;
len_plot = min(5000, len_min);
t_plot = t_comp(1:len_plot);

for k = 1:3
    subplot(3,2,2*k-1);
    plot(t_plot, canales_originales{k}(1:len_plot), 'k', 'LineWidth', 1.5, 'DisplayName', 'Original');
    hold on;
    plot(t_plot, canales_demux{k}(1:len_plot), colores{k}, 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Recuperado');
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    title([nombres_canales{k} ' - Comparación Temporal']);
    legend('show');
    grid on;
    hold off;
    
    % Error temporal
    subplot(3,2,2*k);
    error_temporal = canales_demux{k}(1:len_min) - canales_originales{k}(1:len_min);
    plot(t_comp(1:len_plot), error_temporal(1:len_plot), 'k', 'LineWidth', 1);
    xlabel('Tiempo (s)');
    ylabel('Error');
    title([nombres_canales{k} ' - Error de Recuperación']);
    grid on;
    
    % Estadísticas
    error_rms = sqrt(mean(error_temporal.^2));
    error_max = max(abs(error_temporal));
    text(0.02, 0.95, sprintf('RMS: %.4f\nMáx: %.4f', error_rms, error_max), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', 'BackgroundColor', 'white');
end

%% ANÁLISIS CUANTITATIVO OPTIMIZADO
fprintf('\n=== ANÁLISIS CUANTITATIVO DE RECUPERACIÓN ===\n');

correlaciones = zeros(3,1);
snr_valores = zeros(3,1);
errores_relativos = zeros(3,1);

for k = 1:3
    fprintf('\n%s:\n', nombres_canales{k});
    fprintf('----------------------------------------\n');
    
    % Ajustar a longitud común
    canal_demux_ajustado = canales_demux{k}(1:len_min);
    canal_orig_ajustado = canales_originales{k}(1:len_min);
    
    % Cálculos vectorizados (MÁS EFICIENTE)
    error_temporal = canal_demux_ajustado - canal_orig_ajustado;
    error_rms = sqrt(mean(error_temporal.^2));
    error_max = max(abs(error_temporal));
    
    potencia_orig = sqrt(mean(canal_orig_ajustado.^2));
    potencia_demux = sqrt(mean(canal_demux_ajustado.^2));
    
    error_relativo = error_rms / potencia_orig;
    snr_db = 20*log10(potencia_orig / error_rms);
    
    % Correlación
    R = corrcoef(canal_demux_ajustado, canal_orig_ajustado);
    correlacion_valor = R(1,2);
    
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
    
    % Análisis espectral (calcular FFT solo si no existe)
    X_orig_full = fft(canal_orig_ajustado);
    X_demux_full = fft(canal_demux_ajustado);
    
    energia_orig = sum(abs(X_orig_full).^2);
    energia_demux = sum(abs(X_demux_full).^2);
    
    % Banda de audio
    f_analisis = linspace(0, fs_entrada_orig, len_min);
    idx_audio = (f_analisis >= 300 & f_analisis <= 3400);
    energia_orig_audio = sum(abs(X_orig_full(idx_audio)).^2);
    energia_demux_audio = sum(abs(X_demux_full(idx_audio)).^2);
    
    eficiencia_total = energia_demux / energia_orig;
    eficiencia_audio = energia_demux_audio / energia_orig_audio;
    
    fprintf('  ANÁLISIS ESPECTRAL:\n');
    fprintf('    Eficiencia total: %.4f (%.2f%%)\n', eficiencia_total, eficiencia_total*100);
    fprintf('    Eficiencia banda audio: %.4f (%.2f%%)\n', eficiencia_audio, eficiencia_audio*100);
    
    % Evaluación
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

%% GRÁFICA 4: ANÁLISIS DE ALIASING
figure('Position', [400, 400, 1400, 800]);

for k = 1:3
    subplot(2,3,k);
    
    X_filt = fft(canales_filtrados{k});
    L_filt = length(X_filt);
    magX_filt_dB = 20*log10(abs(X_filt(1:L_filt/2)) + eps);
    f_filt_plot = linspace(0, fs_multiplexada/2, L_filt/2);
    
    plot(f_filt_plot/1000, magX_filt_dB, colores{k}, 'LineWidth', 1.5);
    xlabel('Frecuencia (kHz)');
    ylabel('Magnitud (dB)');
    title([nombres_canales{k} ' - Antes Downsampling']);
    grid on;
    xlim([0, 60]);
    
    hold on;
    ylims_current = ylim;
    for aliasing_freq = fs_entrada_orig/1000:fs_entrada_orig/1000:60
        line([aliasing_freq aliasing_freq], ylims_current, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 0.5);
    end
    
    band_center = (bandas{k}(1) + bandas{k}(2))/2000;
    line([band_center band_center], ylims_current, 'Color', colores{k}, 'LineStyle', '--', 'LineWidth', 2);
    hold off;
    
    % Después downsampling
    subplot(2,3,k+3);
    X_demux = fft(canales_demux{k}(1:len_min));
    X_orig = fft(canales_originales{k}(1:len_min));
    
    magX_demux_dB = 20*log10(abs(X_demux(1:len_min/2)) + eps);
    magX_orig_dB = 20*log10(abs(X_orig(1:len_min/2)) + eps);
    f_comp_plot = linspace(0, fs_entrada_orig/2, len_min/2);
    
    plot(f_comp_plot, magX_orig_dB, colores{k}, 'LineWidth', 2, 'DisplayName', 'Original');
    hold on;
    plot(f_comp_plot, magX_demux_dB, colores{k}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Recuperado');
    
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
fprintf('? Método: Filtrado + Downsampling combinado\n');
fprintf('? Optimizaciones aplicadas:\n');
fprintf('  - Vectorización de operaciones\n');
fprintf('  - Reducción de accesos a memoria\n');
fprintf('  - Cálculos compartidos entre canales\n');
fprintf('  - Pre-asignación de memoria\n');
fprintf('? %d canales procesados exitosamente\n', length(canales_demux));

fprintf('\nESTADÍSTICAS GLOBALES:\n');
fprintf('  Correlación promedio: %.4f\n', mean(correlaciones));
fprintf('  SNR promedio: %.2f dB\n', mean(snr_valores));
fprintf('  Error relativo promedio: %.4f%%\n', mean(errores_relativos)*100);

if mean(correlaciones) > 0.9
    fprintf('  ? RECUPERACIÓN EXITOSA\n');
elseif mean(correlaciones) > 0.7
    fprintf('  ?? RECUPERACIÓN PARCIAL\n');
else
    fprintf('  ? RECUPERACIÓN DEFICIENTE\n');
end

%% GUARDAR RESULTADOS
save('resultados_demultiplexacion_optimizada.mat', 'canales_demux', 'canales_filtrados', 'correlaciones', 'snr_valores', 'errores_relativos');

fprintf('\n? Resultados guardados\n');
fprintf('? Archivos de audio: canal_1_demux.wav, canal_2_demux.wav, canal_3_demux.wav\n');
fprintf('? Análisis completado\n\n');