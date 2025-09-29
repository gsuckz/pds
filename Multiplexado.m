clear; clc;

%% CARGAR SEÑALES
load('archivos_procesados.mat', 'archivos_procesados');

% Parámetros
fs_entrada = 8000;
fs = 120000;
factor = 15;
M = factor;

%% CARGAR CANALES DE AUDIO (UNA SOLA VEZ)
canal1 = audioread(archivos_procesados{1});
canal2 = audioread(archivos_procesados{2});
canal3 = audioread(archivos_procesados{3});

% Precalcular constantes
L = length(canal1);
t = (0:L-1)/fs_entrada;

%% ANÁLISIS ESPECTRAL INDIVIDUAL
canales = {canal1, canal2, canal3};
nombres = {'Canal 1', 'Canal 2', 'Canal 3'};
colores = {'b', 'r', 'g'};

N = L;
f = linspace(0, fs_entrada, N);
f_nyquist = f(1:N/2);

%% GRÁFICA 1: Formas de onda
figure('Position', [100, 100, 1200, 800]);
for k = 1:3
    subplot(3,1,k);
    plot(t, canales{k}, colores{k}, 'LineWidth', 1);
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    title([nombres{k} ' - Forma de onda temporal']);
    grid on;
    
    rms_val = sqrt(mean(canales{k}.^2));
    max_val = max(abs(canales{k}));
    fprintf('%s - RMS: %.4f, Máximo: %.4f\n', nombres{k}, rms_val, max_val);
end

%% GRÁFICA 2: Espectros individuales
figure('Position', [200, 200, 1200, 800]);
for k = 1:3
    X = fft(canales{k});
    magX_dB = 20*log10(abs(X(1:N/2)) + eps);
    
    subplot(3,1,k);
    plot(f_nyquist, magX_dB, colores{k}, 'LineWidth', 1.5);
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title([nombres{k} ' - Espectro de frecuencia (0-4kHz)']);
    grid on;
    xlim([0, 4000]);
    ylim([min(magX_dB)-10, max(magX_dB)+5]);
end

%% GRÁFICA 3: Comparación espectral
figure('Position', [300, 300, 1200, 600]);

subplot(2,1,1);
hold on;
for k = 1:3
    X = fft(canales{k});
    magX_dB = 20*log10(abs(X(1:N/2)) + eps);
    plot(f_nyquist, magX_dB, colores{k}, 'LineWidth', 2, 'DisplayName', nombres{k});
end
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Comparación espectral de las 3 señales de entrada');
legend('show', 'Location', 'best');
grid on;
xlim([0, 4000]);
hold off;

%% ANÁLISIS CON DESPLAZAMIENTO EN FRECUENCIA
subplot(2,1,2);
hold on;

% Vectorización del desplazamiento para canales 1 y 3
mult_desp = (-1).^(0:L-1)';  % Vector de alternancia precalculado

for k = [1, 3]
    canal_desplazado = canales{k} .* mult_desp;  % Vectorizado
    X_desp = fft(canal_desplazado);
    magX_desp_dB = 20*log10(abs(X_desp(1:N/2)) + eps);
    plot(f_nyquist, magX_desp_dB, colores{k}, 'LineWidth', 2, 'LineStyle', '--', ...
         'DisplayName', [nombres{k} ' (desplazado +4kHz)']);
end

% Canal 2 sin desplazamiento
X2 = fft(canales{2});
magX2_dB = 20*log10(abs(X2(1:N/2)) + eps);
plot(f_nyquist, magX2_dB, colores{2}, 'LineWidth', 2, 'DisplayName', [nombres{2} ' (sin desplazamiento)']);

xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Efecto del desplazamiento en frecuencia (+4kHz por alternancia de signo)');
legend('show', 'Location', 'best');
grid on;
xlim([0, 4000]);
hold off;

%% DISEÑO DE FILTROS
banda1 = [12300 15400];
banda2 = [16300 19400];
banda3 = [20300 23400];
bandas = {banda1, banda2, banda3};

ripple_db = 1;
atten_db = 40;

delta_p = (10^(ripple_db/20) - 1)/(10^(ripple_db/20) + 1);
delta_s = 10^(-atten_db/20);

filtro_polifasico = cell(3,1);
h_original = cell(3,1);

for k = 1:3
    f1 = bandas{k}(1)/(fs/2);
    f2 = bandas{k}(2)/(fs/2);
    
    ancho_banda = f2 - f1;
    banda_transicion = 0.1 * ancho_banda;
    
    f_pass = [f1, f2];
    f_stop = [f1 - banda_transicion, f2 + banda_transicion];
    
    [N_filt, Wn, beta, ftype] = kaiserord([f_stop(1), f_pass(1), f_pass(2), f_stop(2)], ...
                                     [0 1 0], [delta_s, delta_p, delta_s]);
    
    N_filt = M * ceil(N_filt/M);
    
    fprintf('Filtro %d: Orden N=%d, beta=%.2f\n', k, N_filt, beta);
    
    h = fir1(N_filt, Wn, ftype, kaiser(N_filt+1, beta), 'noscale');
    h_original{k} = h;
    
    L_h = length(h);
    if mod(L_h, M) ~= 0
        h = [h, zeros(1, M - mod(L_h, M))];
        L_h = length(h);
    end
    
    H_matrix = reshape(h, M, L_h/M);
    filtro_polifasico{k} = H_matrix;
    
    fprintf('Filtro %d descompuesto en %d subfiltros de longitud %d\n', k, M, L_h/M);
end

save('filtro_polifasico.mat', 'filtro_polifasico');

%% MULTIPLEXACIÓN OPTIMIZADA
fprintf('\n=== MULTIPLEXACIÓN POLIFÁSICA OPTIMIZADA ===\n');
fprintf('Procesando %d muestras...\n', L);

% Pre-asignar memoria
L_salida = L * factor;
salidaMux = zeros(L_salida, 1);
salidaSinFiltrar = zeros(L_salida, 1);

% Pre-calcular desplazamiento en frecuencia (VECTORIZADO)
desp_vec = (-1).^(0:L-1)';

% Aplicar desplazamiento a canales 1 y 3
canal1_desp = canal1 .* desp_vec;
canal3_desp = canal3 .* desp_vec;

% Inicializar memorias de filtros
memoria_filtros = cell(3,1);
for k = 1:3
    num_taps = size(filtro_polifasico{k}, 2);
    memoria_filtros{k} = zeros(M, num_taps);
end

tic;
for n = 1:L
    % Progreso cada 20000 muestras (menos I/O)
    if mod(n, 20000) == 0
        fprintf('Procesada muestra %d de %d (%.1f%%)\n', n, L, (n/L)*100);
    end
    
    idx = (n-1)*factor + 1;
    
    % Obtener muestras (ya con desplazamiento aplicado)
    muestras = [canal1_desp(n), canal2(n), canal3_desp(n)];
    
    % Pre-asignar salida filtrada
    salida_filtrada = zeros(factor, 1);
    
    % Procesar cada canal
    for canal_idx = 1:3
        muestra_actual = muestras(canal_idx);
        H_poly = filtro_polifasico{canal_idx};
        memoria = memoria_filtros{canal_idx};
        
        % Procesar subfiltros polifásicos
        for m = 1:M
            % Actualizar memoria (shift optimizado)
            memoria(m, 2:end) = memoria(m, 1:end-1);
            memoria(m, 1) = muestra_actual;
            
            % Producto punto vectorizado
            y_m = sum(H_poly(m, :) .* memoria(m, :));
            salida_filtrada(m) = salida_filtrada(m) + y_m;
        end
        
        memoria_filtros{canal_idx} = memoria;
    end
    
    % Guardar salidas
    salidaMux(idx:idx+factor-1) = 15 * salida_filtrada;
    salidaSinFiltrar(idx:idx+factor-1) = 10 * sum(muestras);
end

salidaMux = 10 * salidaMux;
tiempo_mux = toc;

fprintf('Procesamiento completado en %.3f segundos.\n', tiempo_mux);

save('salidaMux_polifasico.mat','salidaMux');

%% ANÁLISIS DE RESULTADOS
fprintf('\n=== ANÁLISIS ESPECTRAL ===\n');

N_mux = length(salidaMux);
f_mux = linspace(0, fs, N_mux);

% Calcular FFT una sola vez
X_mux = fft(salidaMux);
magX_mux = 20*log10(abs(X_mux));

%% GRÁFICA: Espectro completo
figure;
plot(f_mux, magX_mux);
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Espectro de la señal multiplexada FDM polifásica a 120 kHz');
grid on;

%% GRÁFICA: Hasta Nyquist
figure;
plot(f_mux(1:N_mux/2), magX_mux(1:N_mux/2));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Espectro hasta Nyquist (60 kHz) - Implementación Polifásica');
grid on;

hold on;
ylims = ylim;
for k = 1:3
    line([bandas{k}(1) bandas{k}(1)], ylims, 'Color', 'r', 'LineStyle', '--', 'DisplayName', sprintf('Banda %d inf', k));
    line([bandas{k}(2) bandas{k}(2)], ylims, 'Color', 'r', 'LineStyle', '--', 'DisplayName', sprintf('Banda %d sup', k));
end
hold off;

%% GRÁFICA: Forma de onda
t_mux = (0:N_mux-1)/fs;
figure;
plot(t_mux, salidaMux);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Forma de onda de la señal multiplexada FDM polifásica a 120 kHz');
grid on;

%% ANÁLISIS DE EFICIENCIA
fprintf('\n=== ANÁLISIS DE LA IMPLEMENTACIÓN ===\n');
fprintf('Factor de upsampling: %d\n', factor);
fprintf('Número de canales: 3\n');
fprintf('Longitud de señal de entrada: %d muestras\n', L);
fprintf('Longitud de señal de salida: %d muestras\n', L_salida);
fprintf('Frecuencia de salida: %d Hz\n', fs);
fprintf('Tiempo de procesamiento: %.3f segundos\n', tiempo_mux);

for k = 1:3
    N_filtro = size(filtro_polifasico{k}, 2) * M;
    fprintf('Filtro %d - Orden total: %d, Subfiltros: %dx%d\n', ...
            k, N_filtro, M, size(filtro_polifasico{k}, 2));
end

%% RESPUESTA EN FRECUENCIA DE LOS FILTROS
figure;
for k = 1:3
    [H, w] = freqz(h_original{k}, 1, 1024, fs);
    
    subplot(3,1,k);
    plot(w, 20*log10(abs(H)));
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(sprintf('Respuesta en frecuencia - Filtro Kaiser %d', k));
    grid on;
    
    hold on;
    ylims = ylim;
    line([bandas{k}(1) bandas{k}(1)], ylims, 'Color', 'r', 'LineStyle', '--');
    line([bandas{k}(2) bandas{k}(2)], ylims, 'Color', 'r', 'LineStyle', '--');
    ylim([-60, 5]);
    hold off;
end

%% COMPARACIÓN CON FUNCIÓN FILTER
fprintf('\n=== COMPARACIÓN CON FUNCIÓN FILTER ===\n');
tic;

% Pre-asignar canales con upsampling
canal1_up = zeros(L_salida, 1);
canal2_up = zeros(L_salida, 1);
canal3_up = zeros(L_salida, 1);

% Upsampling optimizado (asignación directa)
idx_up = 1:factor:L_salida;
canal1_up(idx_up) = canal1_desp;
canal2_up(idx_up) = canal2;
canal3_up(idx_up) = canal3_desp;

% Filtrar
canal1_filtrado = filter(h_original{1}, 1, canal1_up);
canal2_filtrado = filter(h_original{2}, 1, canal2_up);
canal3_filtrado = filter(h_original{3}, 1, canal3_up);

% Combinar
salidaMux_filter = 10 * (canal1_filtrado + canal2_filtrado + canal3_filtrado);

tiempo_filter = toc;
fprintf('Tiempo de procesamiento con filter(): %.4f segundos\n', tiempo_filter);

%% ANÁLISIS COMPARATIVO
diferencia = salidaMux_filter - salidaMux;
error_rms = sqrt(mean(diferencia.^2));
error_max = max(abs(diferencia));
error_relativo = error_rms / sqrt(mean(salidaMux.^2));

fprintf('Error RMS entre métodos: %.2e\n', error_rms);
fprintf('Error máximo absoluto: %.2e\n', error_max);
fprintf('Error relativo (%%): %.6f%%\n', error_relativo * 100);

%% COMPARACIÓN ESPECTRAL
X_filter = fft(salidaMux_filter);
X_polifasico = fft(salidaMux);

magX_filter_dB = 20*log10(abs(X_filter));
magX_polifasico_dB = 20*log10(abs(X_polifasico));

%% GRÁFICAS COMPARATIVAS
figure;
subplot(2,2,1);
plot(f_mux(1:N_mux/2)/1000, magX_filter_dB(1:N_mux/2), 'b-', 'LineWidth', 1.5);
hold on;
plot(f_mux(1:N_mux/2)/1000, magX_polifasico_dB(1:N_mux/2), 'r--', 'LineWidth', 1);
xlabel('Frecuencia (kHz)');
ylabel('Magnitud (dB)');
title('Comparación: filter() vs Polifásico');
legend('filter()', 'Polifásico', 'Location', 'best');
grid on;

subplot(2,2,2);
t_comp = (0:min(2000, N_mux-1))/fs;
plot(t_comp, salidaMux_filter(1:length(t_comp)), 'b-', 'LineWidth', 1.5);
hold on;
plot(t_comp, salidaMux(1:length(t_comp)), 'r--', 'LineWidth', 1);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Comparación temporal');
legend('filter()', 'Polifásico');
grid on;

subplot(2,2,3);
plot(t_comp, diferencia(1:length(t_comp)), 'g-', 'LineWidth', 1);
xlabel('Tiempo (s)');
ylabel('Error');
title(sprintf('Error entre métodos (RMS=%.2e)', error_rms));
grid on;

subplot(2,2,4);
diff_espectral = magX_filter_dB(1:N_mux/2) - magX_polifasico_dB(1:N_mux/2);
plot(f_mux(1:N_mux/2)/1000, diff_espectral, 'g-', 'LineWidth', 1);
xlabel('Frecuencia (kHz)');
ylabel('Diferencia (dB)');
title('Diferencia espectral');
grid on;

%% CONCLUSIONES
fprintf('\n=== RESULTADOS COMPARATIVOS ===\n');
if error_relativo < 1e-10
    fprintf('? Los métodos son EQUIVALENTES (error < 1e-8%%)\n');
elseif error_relativo < 1e-6
    fprintf('? Los métodos son PRÁCTICAMENTE IDÉNTICOS (error < 1e-4%%)\n');
else
    fprintf('? Los métodos tienen diferencias (error = %.2e%%)\n', error_relativo*100);
end

fprintf('Mejora de velocidad: %.2fx más rápido que versión original\n', tiempo_filter/tiempo_mux);

save('comparacion_filter.mat', 'salidaMux_filter', 'diferencia', 'error_rms', 'tiempo_filter');

%% ANÁLISIS DE CANALES FILTRADOS INDIVIDUALES
fprintf('\n=== ANÁLISIS DE CANALES FILTRADOS INDIVIDUALMENTE ===\n');

canales_filtrados = cell(3,1);
nombres_canales = {'Canal 1 (12.3-15.4 kHz)', 'Canal 2 (16.3-19.4 kHz)', 'Canal 3 (20.3-23.4 kHz)'};

for canal_idx = 1:3
    fprintf('Procesando %s...\n', nombres_canales{canal_idx});
    
    salida_canal = zeros(L_salida, 1);
    num_taps = size(filtro_polifasico{canal_idx}, 2);
    memoria_canal = zeros(factor, num_taps);
    
    % Seleccionar canal de entrada con desplazamiento ya aplicado
    if canal_idx == 1
        canal_entrada = canal1_desp;
    elseif canal_idx == 2
        canal_entrada = canal2;
    else
        canal_entrada = canal3_desp;
    end
    
    % Procesar
    for n = 1:L
        idx = (n-1)*factor + 1;
        muestra_actual = canal_entrada(n);
        
        salida_filtrada = zeros(factor, 1);
        H_poly = filtro_polifasico{canal_idx};
        
        for m = 1:factor
            memoria_canal(m, 2:end) = memoria_canal(m, 1:end-1);
            memoria_canal(m, 1) = muestra_actual;
            y_m = sum(H_poly(m, :) .* memoria_canal(m, :));
            salida_filtrada(m) = y_m;
        end
        
        salida_canal(idx:idx+factor-1) = salida_filtrada;
    end
    
    canales_filtrados{canal_idx} = salida_canal;
end

% Verificación
suma_verificacion = canales_filtrados{1} + canales_filtrados{2} + canales_filtrados{3};
error_verificacion = max(abs(suma_verificacion - salidaMux));
fprintf('Error de verificación: %.2e\n', error_verificacion);

%% ANÁLISIS ESPECTRAL INDIVIDUAL
fprintf('\n--- Análisis espectral individual ---\n');

N_sep = L_salida;
f_sep = linspace(0, fs, N_sep);
f_plot_sep = f_sep(1:N_sep/2);

espectros = cell(3,1);
bandas_obj = {[12.3 15.4], [16.3 19.4], [20.3 23.4]};

for k = 1:3
    X_canal = fft(canales_filtrados{k});
    espectros{k} = 20*log10(abs(X_canal(1:N_sep/2)) + eps);
    
    banda_obj = [bandas_obj{k}(1)*1000, bandas_obj{k}(2)*1000];
    idx_banda = (f_plot_sep >= banda_obj(1) & f_plot_sep <= banda_obj(2));
    
    potencia_banda = sum(abs(X_canal(idx_banda)).^2);
    potencia_total = sum(abs(X_canal(1:N_sep/2)).^2);
    porcentaje_banda = 100 * potencia_banda / potencia_total;
    
    fprintf('Canal %d: %.1f%% de potencia en banda objetivo\n', k, porcentaje_banda);
end

%% GRÁFICAS INDIVIDUALES
figure('Position', [100, 100, 1400, 900]);

subplot(3,2,1);
hold on;
for k = 1:3
    plot(f_plot_sep/1000, espectros{k}, colores{k}, 'LineWidth', 2, 'DisplayName', nombres_canales{k});
end

for k = 1:3
    ylims = ylim;
    fill([bandas_obj{k}(1) bandas_obj{k}(2) bandas_obj{k}(2) bandas_obj{k}(1)], ...
         [ylims(1) ylims(1) ylims(2) ylims(2)], colores{k}, ...
         'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

xlabel('Frecuencia (kHz)');
ylabel('Magnitud (dB)');
title('Espectros de canales filtrados individualmente');
legend('show', 'Location', 'best');
grid on;
xlim([0, 30]);
hold off;

subplot(3,2,2);
hold on;
for k = 1:3
    f_zoom_start = bandas_obj{k}(1) - 2;
    f_zoom_end = bandas_obj{k}(2) + 2;
    idx_zoom = (f_plot_sep/1000 >= f_zoom_start & f_plot_sep/1000 <= f_zoom_end);
    
    plot(f_plot_sep(idx_zoom)/1000, espectros{k}(idx_zoom), colores{k}, ...
         'LineWidth', 2, 'DisplayName', nombres_canales{k});
    
    line([bandas_obj{k}(1) bandas_obj{k}(1)], ylim, 'Color', colores{k}, 'LineStyle', '--');
    line([bandas_obj{k}(2) bandas_obj{k}(2)], ylim, 'Color', colores{k}, 'LineStyle', '--');
end
xlabel('Frecuencia (kHz)');
ylabel('Magnitud (dB)');
title('Zoom en bandas objetivo');
legend('show', 'Location', 'best');
grid on;
hold off;

for k = 1:3
    subplot(3,2,k+2);
    plot(f_plot_sep/1000, espectros{k}, colores{k}, 'LineWidth', 1.5);
    
    hold on;
    ylims = ylim;
    fill([bandas_obj{k}(1) bandas_obj{k}(2) bandas_obj{k}(2) bandas_obj{k}(1)], ...
         [ylims(1) ylims(1) ylims(2) ylims(2)], colores{k}, ...
         'FaceAlpha', 0.2, 'EdgeColor', colores{k}, 'LineStyle', '--', 'LineWidth', 2);
    
    xlabel('Frecuencia (kHz)');
    ylabel('Magnitud (dB)');
    title(nombres_canales{k});
    grid on;
    xlim([0, 30]);
    
    text(bandas_obj{k}(1) + (bandas_obj{k}(2)-bandas_obj{k}(1))/2, ylims(2)-10, ...
         sprintf('%.1f-%.1f kHz', bandas_obj{k}(1), bandas_obj{k}(2)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', colores{k});
    hold off;
end

%% ANÁLISIS TEMPORAL
figure('Position', [200, 200, 1400, 600]);

t_temp = (0:min(5000, N_sep-1))/fs;
for k = 1:3
    subplot(1,3,k);
    plot(t_temp, canales_filtrados{k}(1:length(t_temp)), colores{k}, 'LineWidth', 1);
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    title([nombres_canales{k} ' - Forma de onda']);
    grid on;
    
    rms_canal = sqrt(mean(canales_filtrados{k}.^2));
    max_canal = max(abs(canales_filtrados{k}));
    text(0.02, 0.95, sprintf('RMS: %.4f\nMáx: %.4f', rms_canal, max_canal), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'white', 'EdgeColor', colores{k});
end

%% ESTADÍSTICAS DETALLADAS
fprintf('\n--- Estadísticas de canales filtrados ---\n');
for k = 1:3
    canal_data = canales_filtrados{k};
    
    fprintf('\n%s:\n', nombres_canales{k});
    fprintf('  RMS: %.6f\n', sqrt(mean(canal_data.^2)));
    fprintf('  Máximo: %.6f\n', max(abs(canal_data)));
    fprintf('  Potencia total: %.6f\n', sum(canal_data.^2)/length(canal_data));
    
    X_canal = fft(canal_data);
    
    banda_target = [bandas_obj{k}(1)*1000, bandas_obj{k}(2)*1000];
    idx_target = (f_plot_sep >= banda_target(1) & f_plot_sep <= banda_target(2));
    potencia_target = sum(abs(X_canal(idx_target)).^2);
    
    idx_fuera = (f_plot_sep < banda_target(1) | f_plot_sep > banda_target(2));
    potencia_fuera = sum(abs(X_canal(idx_fuera)).^2);
    
    relacion_dentro_fuera = 10*log10(potencia_target / potencia_fuera);
    fprintf('  Relación dentro/fuera de banda: %.1f dB\n', relacion_dentro_fuera);
end

save('canales_filtrados_individuales.mat', 'canales_filtrados', 'nombres_canales');
save('filtros_originales.mat', 'h_original');
fprintf('\n? Canales guardados en: canales_filtrados_individuales.mat\n');
fprintf('? Filtros guardados en: filtros_originales.mat\n');