clear; clc;

%% CARGAR DATOS Y CONFIGURACIÓN
fprintf('=== DEMULTIPLEXOR FDM OPTIMIZADO - FILTRADO + DOWNSAMPLING EN UN PASO ===\n\n');

% Cargar la seal multiplexada y datos originales
load('salidaMux_polifasico.mat', 'salidaMux');
load('archivos_procesados.mat', 'archivos_procesados');
load('filtros_originales.mat', 'h_original');

% Cargar seales originales para comparación
canal1_orig = audioread(archivos_procesados{1});
canal2_orig = audioread(archivos_procesados{2});
canal3_orig = audioread(archivos_procesados{3});

% Parámetros del sistema
fs_entrada_orig = 8000;   % Frecuencia original de cada canal
fs_multiplexada = 120000; % Frecuencia de la seal multiplexada
factor = 15;              % Factor de downsampling (120k/8k = 15)

% Definir bandas de frecuencia de cada canal
banda1 = [12300 15400];
banda2 = [16300 19400]; 
banda3 = [20300 23400];
bandas = {banda1, banda2, banda3};
nombres_canales = {'Canal 1', 'Canal 2', 'Canal 3'};
colores = {'b', 'r', 'g'};

% Información de las seales
L_mux = length(salidaMux);
L_orig_esperado = floor(L_mux / factor);

fprintf('Longitud seal multiplexada: %d muestras\n', L_mux);
fprintf('Longitud esperada seales demultiplexadas: %d muestras\n', L_orig_esperado);
fprintf('Factor de downsampling: %d\n\n', factor);

%% CARGAR FILTROS
fprintf('--- Cargando filtros pasabanda ---\n');
h_filtros = h_original;

for canal = 1:3
    fprintf('Canal %d: Filtro de longitud %d coeficientes\n', canal, length(h_filtros{canal}));
end
fprintf('\n');

%% ANÁLISIS ESPECTRAL DE LA SEAL MULTIPLEXADA
fprintf('--- Analizando seal multiplexada de entrada ---\n');

N_mux = length(salidaMux);
f_mux = linspace(0, fs_multiplexada, N_mux);
f_mux_plot = f_mux(1:N_mux/2);

% Calcular espectro de la seal multiplexada
X_mux = fft(salidaMux);
magX_mux_dB = 20*log10(abs(X_mux(1:N_mux/2)) + eps);

% Vector de tiempo para graficar
t_mux = (0:N_mux-1)/fs_multiplexada;

%% GRÁFICA 1: SEAL MULTIPLEXADA DE ENTRADA
figure('Position', [100, 100, 1400, 600]);

% Gráfico en el tiempo
subplot(1,2,1);
muestras_plot = min(5000, length(salidaMux));
plot(t_mux(1:muestras_plot), salidaMux(1:muestras_plot), 'k', 'LineWidth', 1);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Seal Multiplexada FDM - Dominio del Tiempo');
grid on;

% Gráfico en frecuencia
subplot(1,2,2);
plot(f_mux_plot/1000, magX_mux_dB, 'k', 'LineWidth', 1.5);
xlabel('Frecuencia (kHz)');
ylabel('Magnitud (dB)');
title('Espectro de la Seal Multiplexada');
grid on;
xlim([0, 60]);

% Resaltar las bandas de cada canal
hold on;
for canal = 1:3
    limites_y = ylim;
    % Crear rectángulo de color para cada banda
    fill([bandas{canal}(1)/1000, bandas{canal}(2)/1000, bandas{canal}(2)/1000, bandas{canal}(1)/1000], ...
         [limites_y(1), limites_y(1), limites_y(2), limites_y(2)], ...
         colores{canal}, 'FaceAlpha', 0.2, 'EdgeColor', colores{canal}, 'LineStyle', '--');
    
    % Etiqueta del canal
    frecuencia_central = (bandas{canal}(1) + bandas{canal}(2))/2000;
    text(frecuencia_central, limites_y(2)-5, sprintf('C%d', canal), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', colores{canal});
end
hold off;

%% DEMULTIPLEXACIÓN OPTIMIZADA
fprintf('--- Iniciando demultiplexación optimizada ---\n');

% Inicializar contenedores para los resultados
canales_filtrados = cell(3,1);  % Para análisis espectral posterior
canales_demux = cell(3,1);      % Seales demultiplexadas finales

% Procesar cada canal individualmente
for canal_idx = 1:3
    fprintf('Procesando %s (Banda: %.1f - %.1f kHz)...\n', ...
            nombres_canales{canal_idx}, bandas{canal_idx}(1)/1000, bandas{canal_idx}(2)/1000);
    
    % Obtener el filtro correspondiente a este canal
    filtro_actual = h_filtros{canal_idx};
    longitud_filtro = length(filtro_actual);
    longitud_entrada = length(salidaMux);
    
    % Calcular el número de muestras de salida después del downsampling
    longitud_salida = floor(longitud_entrada / factor);
    
    % Inicializar vector de salida
    senal_demux = zeros(longitud_salida, 1);
    
    fprintf('  Aplicando filtrado + downsampling combinado...\n');
    
    % ALGORITMO OPTIMIZADO: Calcular solo las muestras que se van a conservar
    for muestra_salida = 1:longitud_salida
        
        % Calcular la posición en la seal original que corresponde 
        % a esta muestra de salida después del downsampling
        posicion_original = (muestra_salida - 1) * factor + 1;
        
        % Realizar la convolución solo para este punto específico
        suma_convolucion = 0;
        
        for coef_filtro = 1:longitud_filtro
            % Índice en la seal de entrada para este coeficiente del filtro
            indice_entrada = posicion_original - coef_filtro + 1;
            
            % Verificar que el índice esté dentro de los límites válidos
            if indice_entrada >= 1 && indice_entrada <= longitud_entrada
                suma_convolucion = suma_convolucion + filtro_actual(coef_filtro) * salidaMux(indice_entrada);
            end
        end
        
        % Guardar el resultado de la convolución
        senal_demux(muestra_salida) = suma_convolucion;
    end
    
    % APLICAR CORRECCIÓN DE DESPLAZAMIENTO EN FRECUENCIA 
    % (Solo necesaria para canales 1 y 3 debido al proceso de modulación)
    if canal_idx == 1 || canal_idx == 3
        fprintf('  Aplicando corrección de desplazamiento de frecuencia...\n');
        
        % Multiplicar por secuencia alternante (+1, -1, +1, -1, ...)
        multiplicador = 1;
        for muestra = 1:length(senal_demux)
            senal_demux(muestra) = senal_demux(muestra) * multiplicador;
            multiplicador = multiplicador * (-1);  % Alternar signo
        end
    end
    
    % Guardar la seal demultiplexada
    canales_demux{canal_idx} = senal_demux;
    
    % Para análisis posterior: calcular también la seal filtrada completa
    % (esto es solo para las gráficas de análisis espectral)
    senal_filtrada_completa = filter(filtro_actual, 1, salidaMux);
    canales_filtrados{canal_idx} = senal_filtrada_completa;
    
    fprintf('  Canal procesado: %d muestras obtenidas\n', length(senal_demux));
    fprintf('  (Esperado: %d muestras)\n\n', length(canal1_orig));
end

%% GUARDAR ARCHIVOS DE AUDIO
fprintf('--- Guardando archivos de audio recuperados ---\n');

for canal = 1:3
    % Crear nombre del archivo
    nombre_archivo = sprintf('canal_%d_demux.wav', canal);
    
    % Obtener la seal y normalizarla
    senal_a_guardar = canales_demux{canal};
    
    % Normalización: ajustar al rango [-1, 1]
    valor_maximo = max(abs(senal_a_guardar));
    if valor_maximo > 0
        senal_a_guardar = senal_a_guardar / valor_maximo;
    end
    
    % Guardar como archivo WAV
    audiowrite(nombre_archivo, senal_a_guardar, fs_entrada_orig);
    fprintf('✓ %s guardado como %s\n', nombres_canales{canal}, nombre_archivo);
end
fprintf('\n');

%% PREPARAR DATOS PARA ANÁLISIS ESPECTRAL
fprintf('--- Preparando análisis espectral ---\n');

% Parámetros para el análisis de frecuencia
longitud_analisis = length(canales_demux{1});
vector_frecuencias = linspace(0, fs_entrada_orig, longitud_analisis);
frecuencias_plot = vector_frecuencias(1:longitud_analisis/2);
vector_tiempo = (0:longitud_analisis-1)/fs_entrada_orig;

% Seales originales para comparación
canales_originales = {canal1_orig, canal2_orig, canal3_orig};

%% GRÁFICA 2: ANÁLISIS PASO A PASO DEL PROCESAMIENTO
figure('Position', [200, 200, 1600, 1200]);

for canal = 1:3
    
    % COLUMNA 1: Respuesta en frecuencia del filtro
    subplot(3, 4, 4*(canal-1) + 1);
    [respuesta_filtro, frecuencias_respuesta] = freqz(h_filtros{canal}, 1, 2048, fs_multiplexada);
    plot(frecuencias_respuesta/1000, 20*log10(abs(respuesta_filtro)), colores{canal}, 'LineWidth', 2);
    xlabel('Frecuencia (kHz)');
    ylabel('Magnitud (dB)');
    title(sprintf('%s - Respuesta del Filtro', nombres_canales{canal}));
    grid on;
    xlim([0, 30]);
    ylim([-80, 10]);
    
    % Resaltar la banda de paso del filtro
    hold on;
    limites_y = ylim;
    fill([bandas{canal}(1)/1000, bandas{canal}(2)/1000, bandas{canal}(2)/1000, bandas{canal}(1)/1000], ...
         [limites_y(1), limites_y(1), limites_y(2), limites_y(2)], ...
         colores{canal}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold off;
    
    % COLUMNA 2: Espectro después del filtrado (antes del downsampling)
    subplot(3, 4, 4*(canal-1) + 2);
    espectro_filtrado = fft(canales_filtrados{canal});
    magnitud_filtrado_dB = 20*log10(abs(espectro_filtrado(1:length(espectro_filtrado)/2)) + eps);
    frecuencias_filtrado = linspace(0, fs_multiplexada, length(espectro_filtrado));
    frecuencias_filtrado_plot = frecuencias_filtrado(1:length(frecuencias_filtrado)/2);
    
    plot(frecuencias_filtrado_plot/1000, magnitud_filtrado_dB, colores{canal}, 'LineWidth', 1.5);
    xlabel('Frecuencia (kHz)');
    ylabel('Magnitud (dB)');
    title('Después del Filtrado (120 kHz)');
    grid on;
    xlim([0, 30]);
    
    % COLUMNA 3: Espectro después del downsampling
    subplot(3, 4, 4*(canal-1) + 3);
    espectro_demux = fft(canales_demux{canal});
    magnitud_demux_dB = 20*log10(abs(espectro_demux(1:longitud_analisis/2)) + eps);
    
    plot(frecuencias_plot, magnitud_demux_dB, colores{canal}, 'LineWidth', 2);
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title('Después del Downsampling (8 kHz)');
    grid on;
    xlim([0, 4000]);
    
    % COLUMNA 4: Comparación con la seal original
    subplot(3, 4, 4*(canal-1) + 4);
    espectro_original = fft(canales_originales{canal}(1:longitud_analisis));
    magnitud_original_dB = 20*log10(abs(espectro_original(1:longitud_analisis/2)) + eps);
    
    plot(frecuencias_plot, magnitud_original_dB, 'k', 'LineWidth', 2, 'DisplayName', 'Original');
    hold on;
    plot(frecuencias_plot, magnitud_demux_dB, colores{canal}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Recuperado');
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title('Comparación Final');
    legend('Location', 'best');
    grid on;
    xlim([0, 4000]);
    hold off;
end

%% GRÁFICA 3: COMPARACIÓN TEMPORAL DETALLADA
figure('Position', [300, 300, 1400, 900]);

for canal = 1:3
    
    % Ajustar longitudes para comparación justa
    longitud_minima = min(length(canales_demux{canal}), length(canales_originales{canal}));
    tiempo_comparacion = (0:longitud_minima-1)/fs_entrada_orig;
    muestras_temporales = min(5000, longitud_minima);
    tiempo_plot = tiempo_comparacion(1:muestras_temporales);
    
    % COMPARACIÓN TEMPORAL
    subplot(3, 2, 2*(canal-1) + 1);
    plot(tiempo_plot, canales_originales{canal}(1:muestras_temporales), 'k', 'LineWidth', 1.5, 'DisplayName', 'Original');
    hold on;
    plot(tiempo_plot, canales_demux{canal}(1:muestras_temporales), colores{canal}, ...
         'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Recuperado');
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    title([nombres_canales{canal} ' - Comparación Temporal']);
    legend('Location', 'best');
    grid on;
    hold off;
    
    % ANÁLISIS DEL ERROR
    subplot(3, 2, 2*(canal-1) + 2);
    error_temporal = canales_demux{canal}(1:longitud_minima) - canales_originales{canal}(1:longitud_minima);
    plot(tiempo_comparacion(1:muestras_temporales), error_temporal(1:muestras_temporales), 'k', 'LineWidth', 1);
    xlabel('Tiempo (s)');
    ylabel('Error de Amplitud');
    title([nombres_canales{canal} ' - Error de Recuperación']);
    grid on;
    
    % Calcular y mostrar estadísticas del error
    error_rms = sqrt(mean(error_temporal.^2));
    error_maximo = max(abs(error_temporal));
    text(0.02, 0.95, sprintf('Error RMS: %.4f\nError Máximo: %.4f', error_rms, error_maximo), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'black');
end

%% ANÁLISIS CUANTITATIVO COMPLETO
fprintf('=== ANÁLISIS CUANTITATIVO DE LA RECUPERACIÓN ===\n');

% Inicializar vectores para estadísticas globales
correlaciones = zeros(3,1);
valores_snr = zeros(3,1);
errores_relativos = zeros(3,1);

for canal = 1:3
    fprintf('\n%s:\n', nombres_canales{canal});
    fprintf('=====================================\n');
    
    % Preparar seales para comparación
    longitud_minima = min(length(canales_demux{canal}), length(canales_originales{canal}));
    senal_recuperada = canales_demux{canal}(1:longitud_minima);
    senal_original = canales_originales{canal}(1:longitud_minima);
    
    % ANÁLISIS DE ERRORES
    error_temporal = senal_recuperada - senal_original;
    error_rms = sqrt(mean(error_temporal.^2));
    error_maximo = max(abs(error_temporal));
    
    % ANÁLISIS DE POTENCIAS
    potencia_original = sqrt(mean(senal_original.^2));
    potencia_recuperada = sqrt(mean(senal_recuperada.^2));
    
    % MÉTRICAS DE CALIDAD
    error_relativo_porcentual = (error_rms / potencia_original) * 100;
    relacion_seal_ruido_dB = 20*log10(potencia_original / error_rms);
    
    % COEFICIENTE DE CORRELACIÓN
    matriz_correlacion = corrcoef(senal_recuperada, senal_original);
    coeficiente_correlacion = matriz_correlacion(1,2);
    
    % Guardar para estadísticas globales
    correlaciones(canal) = coeficiente_correlacion;
    valores_snr(canal) = relacion_seal_ruido_dB;
    errores_relativos(canal) = error_relativo_porcentual/100;
    
    % MOSTRAR RESULTADOS
    fprintf('ANÁLISIS TEMPORAL:\n');
    fprintf('  Error RMS:              %.6f\n', error_rms);
    fprintf('  Error máximo:           %.6f\n', error_maximo);
    fprintf('  Error relativo:         %.2f%%\n', error_relativo_porcentual);
    
    fprintf('\nMÉTRICAS DE FIDELIDAD:\n');
    fprintf('  Correlación:            %.6f\n', coeficiente_correlacion);
    fprintf('  SNR:                    %.2f dB\n', relacion_seal_ruido_dB);
    fprintf('  Potencia original:      %.6f\n', potencia_original);
    fprintf('  Potencia recuperada:    %.6f\n', potencia_recuperada);
    fprintf('  Relación de potencias:  %.4f\n', potencia_recuperada/potencia_original);
    
    % ANÁLISIS ESPECTRAL
    espectro_original = fft(senal_original);
    espectro_recuperado = fft(senal_recuperada);
    
    energia_total_original = sum(abs(espectro_original).^2);
    energia_total_recuperada = sum(abs(espectro_recuperado).^2);
    
    % Energía en la banda de audio útil (300-3400 Hz)
    frecuencias_analisis = linspace(0, fs_entrada_orig, longitud_minima);
    indices_banda_audio = find(frecuencias_analisis >= 300 & frecuencias_analisis <= 3400);
    energia_audio_original = sum(abs(espectro_original(indices_banda_audio)).^2);
    energia_audio_recuperada = sum(abs(espectro_recuperado(indices_banda_audio)).^2);
    
    eficiencia_total = energia_total_recuperada / energia_total_original;
    eficiencia_banda_audio = energia_audio_recuperada / energia_audio_original;
    
    fprintf('\nANÁLISIS ESPECTRAL:\n');
    fprintf('  Eficiencia energética total:    %.2f%% \n', eficiencia_total*100);
    fprintf('  Eficiencia en banda de audio:   %.2f%% \n', eficiencia_banda_audio*100);
    
    % EVALUACIÓN CUALITATIVA
    if coeficiente_correlacion > 0.95 && error_relativo_porcentual < 5
        evaluacion = 'EXCELENTE';
    elseif coeficiente_correlacion > 0.90 && error_relativo_porcentual < 10
        evaluacion = 'BUENA';
    elseif coeficiente_correlacion > 0.80 && error_relativo_porcentual < 20
        evaluacion = 'ACEPTABLE';
    else
        evaluacion = 'DEFICIENTE';
    end
    
    fprintf('\nEVALUACIÓN GLOBAL:        %s\n', evaluacion);
end

%% GRÁFICA 4: ANÁLISIS DE ALIASING Y COMPONENTES ESPECTRALES
figure('Position', [400, 400, 1400, 800]);

for canal = 1:3
    
    % PARTE SUPERIOR: Espectro antes del downsampling
    subplot(2, 3, canal);
    
    espectro_antes_down = fft(canales_filtrados{canal});
    magnitud_antes_down_dB = 20*log10(abs(espectro_antes_down(1:length(espectro_antes_down)/2)) + eps);
    frecuencias_antes_down = linspace(0, fs_multiplexada, length(espectro_antes_down));
    frecuencias_antes_down_plot = frecuencias_antes_down(1:length(frecuencias_antes_down)/2);
    
    plot(frecuencias_antes_down_plot/1000, magnitud_antes_down_dB, colores{canal}, 'LineWidth', 1.5);
    xlabel('Frecuencia (kHz)');
    ylabel('Magnitud (dB)');
    title([nombres_canales{canal} ' - Espectro Pre-Downsampling']);
    grid on;
    xlim([0, 60]);
    
    % Marcar múltiplos de la frecuencia de muestreo final (posible aliasing)
    hold on;
    for multiplo_fs = fs_entrada_orig/1000:fs_entrada_orig/1000:60
        limites_y = ylim;
        line([multiplo_fs multiplo_fs], limites_y, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 0.5);
    end
    
    % Marcar frecuencia central de la banda original
    frecuencia_central = (bandas{canal}(1) + bandas{canal}(2))/2000;
    limites_y = ylim;
    line([frecuencia_central frecuencia_central], limites_y, 'Color', colores{canal}, 'LineStyle', '--', 'LineWidth', 2);
    hold off;
    
    % PARTE INFERIOR: Comparación final con análisis de diferencias
    subplot(2, 3, canal + 3);
    longitud_minima = min(length(canales_demux{canal}), length(canales_originales{canal}));
    espectro_demux = fft(canales_demux{canal}(1:longitud_minima));
    espectro_orig = fft(canales_originales{canal}(1:longitud_minima));
    
    magnitud_demux_dB = 20*log10(abs(espectro_demux(1:longitud_minima/2)) + eps);
    magnitud_orig_dB = 20*log10(abs(espectro_orig(1:longitud_minima/2)) + eps);
    frecuencias_comp = linspace(0, fs_entrada_orig, longitud_minima);
    frecuencias_comp_plot = frecuencias_comp(1:longitud_minima/2);
    
    plot(frecuencias_comp_plot, magnitud_orig_dB, colores{canal}, 'LineWidth', 2, 'DisplayName', 'Original');
    hold on;
    plot(frecuencias_comp_plot, magnitud_demux_dB, colores{canal}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Recuperado');
    
    % Mostrar diferencia espectral (desplazada para visualización)
    diferencia_espectral_dB = magnitud_demux_dB - magnitud_orig_dB;
    plot(frecuencias_comp_plot, diferencia_espectral_dB - 20, 'r:', 'LineWidth', 1, 'DisplayName', 'Diferencia (dB-20)');
    
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title([nombres_canales{canal} ' - Comparación Espectral Final']);
    legend('Location', 'best');
    grid on;
    xlim([0, 4000]);
    hold off;
end

%% RESUMEN EJECUTIVO
fprintf('\n=== RESUMEN EJECUTIVO DEL DEMULTIPLEXOR OPTIMIZADO ===\n');
fprintf('✓ Método implementado: Filtrado + Downsampling combinado\n');
fprintf('✓ Ventajas del algoritmo optimizado:\n');
fprintf('  • Menor consumo de memoria RAM\n');
fprintf('  • Mayor eficiencia computacional\n');
fprintf('  • Cálculo selectivo (solo muestras necesarias)\n');
fprintf('  • Resultados equivalentes al método convencional\n');
fprintf('✓ Canales procesados exitosamente: %d\n', length(canales_demux));

% ESTADÍSTICAS FINALES CONSOLIDADAS
fprintf('\n--- ESTADÍSTICAS CONSOLIDADAS ---\n');
fprintf('Correlación promedio:     %.4f\n', mean(correlaciones));
fprintf('SNR promedio:             %.2f dB\n', mean(valores_snr));
fprintf('Error relativo promedio:  %.2f%%\n', mean(errores_relativos)*100);

% EVALUACIÓN FINAL DEL SISTEMA
correlacion_promedio = mean(correlaciones);
if correlacion_promedio > 0.9
    resultado_final = '✅ RECUPERACIÓN EXITOSA';
    descripcion = 'Las seales se han recuperado con alta fidelidad';
elseif correlacion_promedio > 0.7
    resultado_final = '⚠️  RECUPERACIÓN PARCIAL';
    descripcion = 'Las seales presentan cierta degradación aceptable';
else
    resultado_final = '❌ RECUPERACIÓN DEFICIENTE';
    descripcion = 'Las seales muestran distorsión significativa';
end

fprintf('\n%s\n', resultado_final);
fprintf('%s\n', descripcion);

%% GUARDAR RESULTADOS FINALES
save('resultados_demultiplexacion_optimizada.mat', ...
     'canales_demux', 'canales_filtrados', 'correlaciones', 'valores_snr', 'errores_relativos');

fprintf('\n✓ Datos guardados en: resultados_demultiplexacion_optimizada.mat\n');
fprintf('✓ Archivos de audio: canal_1_demux.wav, canal_2_demux.wav, canal_3_demux.wav\n');
fprintf('✓ Procesamiento completado exitosamente\n\n');