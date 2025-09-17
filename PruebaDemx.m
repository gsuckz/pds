clear; clc;

%% CARGAR DATOS Y CONFIGURACI�N
fprintf('=== DEMULTIPLEXOR FDM CON FUNCI�N FILTER ===\n\n');

% Cargar la se�al multiplexada y datos originales
load('salidaMux_polifasico.mat', 'salidaMux');
load('archivos_procesados.mat', 'archivos_procesados');
load('filtros_originales.mat', 'h_original');   % << Filtros originales

% Cargar se�ales originales para comparaci�n
canal1_orig = audioread(archivos_procesados{1});
canal2_orig = audioread(archivos_procesados{2});
canal3_orig = audioread(archivos_procesados{3});

% Par�metros
fs_entrada_orig = 8000;   % Frecuencia original
fs_multiplexada = 120000; % Frecuencia de la se�al multiplexada
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

fprintf('Longitud se�al multiplexada: %d muestras\n', L_mux);
fprintf('Longitud esperada se�ales originales: %d muestras\n', L_orig);
fprintf('Factor de downsampling: %d\n\n', factor);

%% USAR FILTROS ORIGINALES
fprintf('--- Usando filtros originales guardados ---\n');
h_filtros = h_original;

for k = 1:3
    fprintf('Filtro %d cargado - Longitud: %d\n', k, length(h_filtros{k}));
end

%% AN�LISIS ESPECTRAL DE LA SE�AL MULTIPLEXADA
fprintf('\n--- An�lisis de la se�al multiplexada ---\n');

N_mux = length(salidaMux);
f_mux = linspace(0, fs_multiplexada, N_mux);
f_mux_plot = f_mux(1:N_mux/2);

X_mux = fft(salidaMux);
magX_mux_dB = 20*log10(abs(X_mux(1:N_mux/2)) + eps);

% Crear vector de tiempo para la se�al multiplexada
t_mux = (0:N_mux-1)/fs_multiplexada;

%% GR�FICA 1: SE�AL MULTIPLEXADA DE ENTRADA
figure('Position', [100, 100, 1400, 600]);

subplot(1,2,1);
plot(t_mux(1:min(5000, end)), salidaMux(1:min(5000, end)), 'k', 'LineWidth', 1);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Se�al Multiplexada FDM - Entrada al Demultiplexor');
grid on;

subplot(1,2,2);
plot(f_mux_plot/1000, magX_mux_dB, 'k', 'LineWidth', 1.5);
xlabel('Frecuencia (kHz)');
ylabel('Magnitud (dB)');
title('Espectro de la Se�al Multiplexada (120 kHz)');
grid on;
xlim([0, 60]);

% Resaltar las bandas de inter�s
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

%% DEMULTIPLEXACI�N USANDO FUNCI�N FILTER
fprintf('\n--- Iniciando demultiplexaci�n con funci�n filter() ---\n');

% Inicializar se�ales demultiplexadas
canales_filtrados = cell(3,1);
canales_demux = cell(3,1);

% Procesar cada canal
for canal_idx = 1:3
    fprintf('Procesando %s (banda: %.1f-%.1f kHz)...\n', ...
            nombres_canales{canal_idx}, bandas{canal_idx}(1)/1000, bandas{canal_idx}(2)/1000);    
    % PASO 1: Filtrar la se�al multiplexada usando filter()
    fprintf('  Aplicando filtro pasabanda...\n');
    senal_filtrada = filter(h_filtros{canal_idx}, 1, salidaMux);
    canales_filtrados{canal_idx} = senal_filtrada;    
    % PASO 2: Downsampling
    fprintf('  Realizando downsampling (factor %d)...\n', factor);
    canal_recuperado = senal_filtrada(1:factor:end);
    
    % PASO 3: Aplicar desplazamiento en frecuencia inverso para canales 1 y 3
    if canal_idx == 1 || canal_idx == 3
        fprintf('  Aplicando correcci�n de desplazamiento en frecuencia...\n');
       multDesplazamiento = 1;
      for n = 1:length(canal_recuperado)
         canal_recuperado(n) = canal_recuperado(n) * multDesplazamiento;
        multDesplazamiento = multDesplazamiento * (-1);
    end
    end
    
    canales_demux{canal_idx} = canal_recuperado;
    
    fprintf('  Longitud canal recuperado: %d muestras\n', length(canal_recuperado));
    fprintf('  Longitud esperada: %d muestras\n\n', length(canal1_orig));
end


for k = 1:3
    % Ajustar nombre de archivo
    nombre_archivo = sprintf('canal%d_recuperado.wav', k);
    
    % Asegurar que los datos estén en rango [-1, 1]
    canal = canales_demux{k};
    canal = canal / max(abs(canal));   % normalización
    
    % Guardar audio
    audiowrite(nombre_archivo, canal, fs_entrada_orig);
    fprintf('Canal %d guardado en %s\n', k, nombre_archivo);
end


%% AN�LISIS ESPECTRAL - PASO A PASO
fprintf('--- An�lisis espectral paso a paso ---\n');

% Par�metros para FFT
N_orig = length(canales_demux{1});
f_orig = linspace(0, fs_entrada_orig, N_orig);
f_orig_plot = f_orig(1:N_orig/2);
t_orig = (0:N_orig-1)/fs_entrada_orig;

%% GR�FICA 2: PROCESO PASO A PASO PARA CADA CANAL
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
    
    % PASO 2: Se�al filtrada (120 kHz)
    subplot(3,4,4*k-2);
    X_filt = fft(canales_filtrados{k});
    magX_filt_dB = 20*log10(abs(X_filt(1:length(X_filt)/2)) + eps);
    f_filt = linspace(0, fs_multiplexada, length(X_filt));
    f_filt_plot = f_filt(1:length(f_filt)/2);
    
    plot(f_filt_plot/1000, magX_filt_dB, colores{k}, 'LineWidth', 1.5);
    xlabel('Frecuencia (kHz)');
    ylabel('Magnitud (dB)');
    title('Despu�s Filtrado (120 kHz)');
    grid on;
    xlim([0, 30]);
    
    % PASO 3: Despu�s del downsampling
    subplot(3,4,4*k-1);
    X_demux = fft(canales_demux{k});
    magX_demux_dB = 20*log10(abs(X_demux(1:N_orig/2)) + eps);
    
    plot(f_orig_plot, magX_demux_dB, colores{k}, 'LineWidth', 2);
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title('Despu�s Downsampling (8 kHz)');
    grid on;
    xlim([0, 4000]);
    
    % PASO 4: Comparaci�n con original
    subplot(3,4,4*k);
    X_orig = fft(canales_originales{k}(1:N_orig));
    magX_orig_dB = 20*log10(abs(X_orig(1:N_orig/2)) + eps);
    
    plot(f_orig_plot, magX_orig_dB, 'k', 'LineWidth', 2, 'DisplayName', 'Original');
    hold on;
    plot(f_orig_plot, magX_demux_dB, colores{k}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Recuperado');
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title('Comparaci�n Final');
    legend('show');
    grid on;
    xlim([0, 4000]);
    hold off;
end




%% GR�FICA 3: COMPARACI�N TEMPORAL
figure('Position', [300, 300, 1400, 900]);

for k = 1:3
    % Ajustar longitudes para comparaci�n
    len_min = min(length(canales_demux{k}), length(canales_originales{k}));
    t_comp = (0:len_min-1)/fs_entrada_orig;
    t_plot = t_comp(1:min(5000, len_min));
    
    subplot(3,2,2*k-1);
    plot(t_plot, canales_originales{k}(1:length(t_plot)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Original');
    hold on;
    plot(t_plot, canales_demux{k}(1:length(t_plot)), colores{k}, 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Recuperado');
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    title([nombres_canales{k} ' - Comparaci�n Temporal']);
    legend('show');
    grid on;
    hold off;
    
    % Error temporal
    subplot(3,2,2*k);
    error_temporal = canales_demux{k}(1:len_min) - canales_originales{k}(1:len_min);
    plot(t_comp(1:min(5000, len_min)), error_temporal(1:min(5000, len_min)), 'k', 'LineWidth', 1);
    xlabel('Tiempo (s)');
    ylabel('Error');
    title([nombres_canales{k} ' - Error de Recuperaci�n']);
    grid on;
    
    % Estad�sticas en el gr�fico
    error_rms = sqrt(mean(error_temporal.^2));
    error_max = max(abs(error_temporal));
    text(0.02, 0.95, sprintf('RMS: %.4f\nM�x: %.4f', error_rms, error_max), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', 'BackgroundColor', 'white');
end

%% AN�LISIS CUANTITATIVO DETALLADO
fprintf('\n=== AN�LISIS CUANTITATIVO DE RECUPERACI�N ===\n');

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
    
    % Correlaci�n
    correlacion = corrcoef(canal_demux_ajustado, canal_orig_ajustado);
    correlacion_valor = correlacion(1,2);
    
    fprintf('  ERRORES TEMPORALES:\n');
    fprintf('    Error RMS: %.6f\n', error_rms);
    fprintf('    Error m�ximo: %.6f\n', error_max);
    fprintf('    Error relativo: %.4f%%\n', error_relativo * 100);
    
    fprintf('  M�TRICAS DE CALIDAD:\n');
    fprintf('    Correlaci�n: %.6f\n', correlacion_valor);
    fprintf('    SNR: %.2f dB\n', snr_db);
    fprintf('    Potencia original: %.6f\n', potencia_orig);
    fprintf('    Potencia recuperada: %.6f\n', potencia_demux);
    fprintf('    Relaci�n potencias: %.4f\n', potencia_demux/potencia_orig);
    
    % An�lisis espectral
    X_orig_full = fft(canal_orig_ajustado);
    X_demux_full = fft(canal_demux_ajustado);
    
    % Energ�a total
    energia_orig = sum(abs(X_orig_full).^2);
    energia_demux = sum(abs(X_demux_full).^2);
    
    % Energ�a en banda de audio (300-3400 Hz)
    f_analisis = linspace(0, fs_entrada_orig, len_min);
    idx_audio = find(f_analisis >= 300 & f_analisis <= 3400);
    energia_orig_audio = sum(abs(X_orig_full(idx_audio)).^2);
    energia_demux_audio = sum(abs(X_demux_full(idx_audio)).^2);
    
    eficiencia_total = energia_demux / energia_orig;
    eficiencia_audio = energia_demux_audio / energia_orig_audio;
    
    fprintf('  AN�LISIS ESPECTRAL:\n');
    fprintf('    Eficiencia total: %.4f (%.2f%%)\n', eficiencia_total, eficiencia_total*100);
    fprintf('    Eficiencia banda audio: %.4f (%.2f%%)\n', eficiencia_audio, eficiencia_audio*100);
    
    % Distorsi�n arm�nica (aproximada)
    % Comparar componentes espectrales principales
    [pks_orig, locs_orig] = findpeaks(abs(X_orig_full(1:len_min/2)), 'MinPeakHeight', max(abs(X_orig_full))*0.1);
    [pks_demux, locs_demux] = findpeaks(abs(X_demux_full(1:len_min/2)), 'MinPeakHeight', max(abs(X_demux_full))*0.1);
    
    fprintf('    Picos espectrales orig: %d\n', length(pks_orig));
    fprintf('    Picos espectrales demux: %d\n', length(pks_demux));
    
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
    
    fprintf('  EVALUACI�N GENERAL: %s\n', calidad);
end

%% GR�FICA 4: AN�LISIS DE ALIASING Y EFECTOS
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
        line([aliasing_freq aliasing_freq], ylim, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 0.5);

    end
    
    % Marcar banda original
    band_center = (bandas{k}(1) + bandas{k}(2))/2000;
    line([band_center band_center], ylim, 'Color', colores{k}, 'LineStyle', '--', 'LineWidth', 2);
    hold off;
    
    % An�lisis despu�s del downsampling
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
    title([nombres_canales{k} ' - Comparaci�n Final']);
    legend('show');
    grid on;
    xlim([0, 4000]);
    hold off;
end

%% RESUMEN FINAL
fprintf('\n=== RESUMEN DE LA DEMULTIPLEXACI�N ===\n');
fprintf('? M�todo: Filtros pasabanda + downsampling directo\n');
fprintf('? Implementaci�n: Funci�n filter() de MATLAB\n');
fprintf('? %d canales procesados\n', length(canales_demux));

% Estad�sticas globales
correlaciones = zeros(3,1);
snr_valores = zeros(3,1);
errores_relativos = zeros(3,1);

for k = 1:3
    len_min = min(length(canales_demux{k}), length(canales_originales{k}));
    canal_demux_ajustado = canales_demux{k}(1:len_min);
    canal_orig_ajustado = canales_originales{k}(1:len_min);
    
    error_temporal = canal_demux_ajustado - canal_orig_ajustado;
    error_rms = sqrt(mean(error_temporal.^2));
    potencia_orig = sqrt(mean(canal_orig_ajustado.^2));
    
    correlacion = corrcoef(canal_demux_ajustado, canal_orig_ajustado);
    correlaciones(k) = correlacion(1,2);
    snr_valores(k) = 20*log10(potencia_orig / error_rms);
    errores_relativos(k) = error_rms / potencia_orig;
end

fprintf('\nESTAD�STICAS GLOBALES:\n');
fprintf('  Correlaci�n promedio: %.4f\n', mean(correlaciones));
fprintf('  SNR promedio: %.2f dB\n', mean(snr_valores));
fprintf('  Error relativo promedio: %.4f%%\n', mean(errores_relativos)*100);

if mean(correlaciones) > 0.9
    fprintf('  ? RECUPERACI�N EXITOSA - Las se�ales se pueden recuperar satisfactoriamente\n');
elseif mean(correlaciones) > 0.7
    fprintf('  ??  RECUPERACI�N PARCIAL - Las se�ales se recuperan con cierta degradaci�n\n');
else
    fprintf('  ? RECUPERACI�N DEFICIENTE - Las se�ales presentan alta distorsi�n\n');
end

%% GUARDAR RESULTADOS
save('resultados_demultiplexacion.mat', 'canales_demux', 'canales_filtrados', 'correlaciones', 'snr_valores', 'errores_relativos');

fprintf('\n? Resultados guardados en: resultados_demultiplexacion.mat\n');
fprintf('? An�lisis completado\n\n');