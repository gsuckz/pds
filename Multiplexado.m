clear; clc;
% Cargar seales
load('archivos_procesados.mat', 'archivos_procesados');
% Parmetros
fs_entrada = 8000;  % Frecuencia de muestreo original de las seales
% Cargar canales de audio
canal1 = audioread(archivos_procesados{1});
canal2 = audioread(archivos_procesados{2});
canal3 = audioread(archivos_procesados{3});
% Crear vector de tiempo
t = (0:length(canal1)-1)/fs_entrada;
%% Anlisis espectral individual
canales = {canal1, canal2, canal3};
nombres = {'Canal 1', 'Canal 2', 'Canal 3'};
% FFT de cada canal
N = length(canal1);
f = linspace(0, fs_entrada, N);
f_nyquist = f(1:N/2);
% Colores para las grficas
colores = {'b', 'r', 'g'};
%% Grfica 1: Formas de onda en el tiempo
figure('Position', [100, 100, 1200, 800]);
for k = 1:3
    subplot(3,1,k);
    plot(t, canales{k}, colores{k}, 'LineWidth', 1);
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    title([nombres{k} ' - Forma de onda temporal']);
    grid on;
    
    % Estadsticas bsicas
    rms_val = sqrt(mean(canales{k}.^2));
    max_val = max(abs(canales{k}));
    fprintf('%s - RMS: %.4f, Mximo: %.4f\n', nombres{k}, rms_val, max_val);
end
%% Grfica 2: Espectros individuales
figure('Position', [200, 200, 1200, 800]);
for k = 1:3
    % Calcular FFT
    X = fft(canales{k});
    magX_dB = 20*log10(abs(X(1:N/2)) + eps);  % +eps para evitar log(0)
    
    subplot(3,1,k);
    plot(f_nyquist, magX_dB, colores{k}, 'LineWidth', 1.5);
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title([nombres{k} ' - Espectro de frecuencia (0-4kHz)']);
    grid on;
    xlim([0, 4000]);  % Mostrar hasta Nyquist
    ylim([min(magX_dB)-10, max(magX_dB)+5]);
    
end
%% Grfica 3: Comparacin espectral superpuesta
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
title('Comparacin espectral de las 3 seales de entrada');
legend('show', 'Location', 'best');
grid on;
xlim([0, 4000]);
hold off;
%% Anlisis despus del desplazamiento en frecuencia
subplot(2,1,2);
hold on;
% Simular el efecto del desplazamiento por alternancia de signo
multDesplazamintoFrecuencia = 1;
for k = [1, 3]  % Solo canales 1 y 3 tienen desplazamiento
    % Crear seal con desplazamiento
    canal_desplazado = zeros(size(canales{k}));
    for n = 1:length(canales{k})
        canal_desplazado(n) = canales{k}(n) * multDesplazamintoFrecuencia;
        multDesplazamintoFrecuencia = multDesplazamintoFrecuencia * (-1);
    end 
    % FFT de la seal desplazada
    X_desp = fft(canal_desplazado);
    magX_desp_dB = 20*log10(abs(X_desp(1:N/2)) + eps);
    plot(f_nyquist, magX_desp_dB, colores{k}, 'LineWidth', 2, 'LineStyle', '--', ...
         'DisplayName', [nombres{k} ' (desplazado +4kHz)']);
    % Reset para siguiente canal
    multDesplazamintoFrecuencia = 1;
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

%% Bandas objetivo del multiplexado
banda1 = [12300 15400];
banda2 = [16300 19400];
banda3 = [20300 23400];
bandas = {banda1, banda2, banda3};

load('archivos_procesados.mat', 'archivos_procesados');
factor = 15;        % Upsampling corregido a 15
M = factor;         % Nmero de subfiltros polifsicos

%Filtros Multiplexado
fs = 120000;
banda1 = [12300 15400];       % Banda de paso
banda2 = [16300 19400];       % Banda de paso  
banda3 = [20300 23400];       % Banda de paso
bandas = {banda1, banda2, banda3};

% Especificaciones para Kaiser
ripple_db = 1;          % 1% = 20*log10(1.01) ? 0.09 dB ripple en banda de paso
atten_db = 40;          % 20 dB atenuacin en banda de rechazo

% Calcular tolerancias
delta_p = (10^(ripple_db/20) - 1)/(10^(ripple_db/20) + 1);  % Tolerancia banda paso
delta_s = 10^(-atten_db/20);                                % Tolerancia banda rechazo

% Disear filtros Kaiser y descomponerlos en subfiltros polifsicos
filtro_polifasico = cell(3,1);


for k = 1:3
    % Normalizar frecuencias a Nyquist
    f1 = bandas{k}(1)/(fs/2);
    f2 = bandas{k}(2)/(fs/2);    
    % Calcular frecuencias de transicin (agregar banda de transicin)
    ancho_banda = f2 - f1;
    banda_transicion = 0.1 * ancho_banda;  % 10% del ancho de banda
    
    f_pass = [f1, f2];
    f_stop = [f1 - banda_transicion, f2 + banda_transicion];
    
    % Usar kaiserord para calcular orden y beta
    [N, Wn, beta, ftype] = kaiserord([f_stop(1), f_pass(1), f_pass(2), f_stop(2)], ...
                                     [0 1 0], [delta_s, delta_p, delta_s]);
    
    % Asegurar que N sea mltiplo de M para la descomposicin polifsica
    N = M * ceil(N/M);
    
    fprintf('Filtro %d: Orden N=%d, beta=%.2f\n', k, N, beta);
    
    % Disear filtro Kaiser
    h = fir1(N, Wn, ftype, kaiser(N+1, beta), 'noscale');
    h_original{k} = h;
    
    % Descomponer en M subfiltros polifsicos
    L = length(h);
    if mod(L, M) ~= 0
        h = [h, zeros(1, M - mod(L, M))];
        L = length(h);
    end
    
    % Matriz polifsica: cada fila es un subfiltro
    H_matrix = reshape(h, M, L/M);
    
    % Guardar subfiltros como vectores fila
    filtro_polifasico{k} = H_matrix;
    
    % Mostrar informacin del filtro
    fprintf('Filtro %d descompuesto en %d subfiltros de longitud %d\n', k, M, L/M);
end

save('filtro_polifasico.mat', 'filtro_polifasico');

% Cargar canales de audio
canal1 = audioread(archivos_procesados{1});
canal2 = audioread(archivos_procesados{2});
canal3 = audioread(archivos_procesados{3});

L = length(canal1); % Todos son iguales

% Inicializar salida y memorias de los filtros
salidaMux = zeros(L*factor, 1);
salidaSinFiltrar = zeros(L*factor, 1);

% Memorias para cada subfiltro de cada canal (estados internos)
memoria_filtros = cell(3,1);
for k = 1:3
    num_taps = size(filtro_polifasico{k}, 2);
    memoria_filtros{k} = zeros(M, num_taps);  % M subfiltros, cada uno con num_taps estados
end

multDesplazamintoFrecuencia = 1;

fprintf('Procesando %d muestras...\n', L);

for n = 1:L
    % Mostrar progreso cada 10000 muestras
    if mod(n, 10000) == 0
        fprintf('Procesada muestra %d de %d\n', n, L);
    end
    
    % Posicin inicial en la seal de salida
    idx = (n-1)*factor + 1;
    
    % Obtener muestras de entrada
    muestra1 = canal1(n) * multDesplazamintoFrecuencia;
    muestra2 = canal2(n);
    muestra3 = canal3(n) * multDesplazamintoFrecuencia;
    multDesplazamintoFrecuencia = multDesplazamintoFrecuencia * (-1);
    
    muestras = [muestra1, muestra2, muestra3];
    
    % Inicializar salidas del filtrado
    salida_filtrada = zeros(factor, 1);
    
    % Procesar cada canal con su filtro polifsico
    for canal = 1:3
        muestra_actual = muestras(canal);
        H_poly = filtro_polifasico{canal};
        memoria = memoria_filtros{canal};
        
        % Procesar cada subfiltro polifsico
        for m = 1:M
            % Desplazar memoria del subfiltro m
            memoria(m, 2:end) = memoria(m, 1:end-1);
            memoria(m, 1) = muestra_actual;
            
            % Calcular salida del subfiltro m
            y_m = sum(H_poly(m, :) .* memoria(m, :));
            
            % Acumular en la posicin correcta de la salida
            salida_filtrada(m) = salida_filtrada(m) + y_m;
        end
        
        % Actualizar memoria
        memoria_filtros{canal} = memoria;
    end
    
    % Guardar las M muestras de salida
    salidaMux(idx:idx+factor-1) = 15*salida_filtrada;
    
    % Salida sin filtrar para comparacin
    suma_sin_filtrar = sum(muestras);
    salidaSinFiltrar(idx:idx+factor-1) = 10*suma_sin_filtrar;
end
salidaMux = 10*salidaMux;
fprintf('Procesamiento completado.\n');

save('salidaMux_polifasico.mat','salidaMux');

%% Anlisis de resultados
fs_salida = 120000;  % frecuencia de muestreo
x = salidaMux;

%% FFT
N = length(x);
X = fft(x);
f = linspace(0, fs_salida, N);

%% Magnitud en dB
magX = 20*log10(abs(X));

%% Graficar espectro completo
figure;
plot(f, magX);
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Espectro de la seal multiplexada FDM polifsica a 120 kHz');
grid on;

%% Graficar solo hasta Nyquist (60 kHz)
figure;
plot(f(1:N/2), magX(1:N/2));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Espectro hasta Nyquist (60 kHz) - Implementacin Polifsica');
grid on;

% Resaltar las bandas de inters
hold on;
ylims = ylim;
for k = 1:3
    % Lneas verticales para R2016a
    line([bandas{k}(1) bandas{k}(1)], ylims, 'Color', 'r', 'LineStyle', '--', 'DisplayName', sprintf('Banda %d inf', k));
    line([bandas{k}(2) bandas{k}(2)], ylims, 'Color', 'r', 'LineStyle', '--', 'DisplayName', sprintf('Banda %d sup', k));
end
hold off;

%% Graficar forma de onda seal multiplexada 
t = (0:N-1)/fs_salida;
figure;
plot(t, x);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Forma de onda de la seal multiplexada FDM polifsica a 120 kHz');
grid on;

%% Anlisis de eficiencia
fprintf('\n=== ANLISIS DE LA IMPLEMENTACIN ===\n');
fprintf('Factor de upsampling: %d\n', factor);
fprintf('Nmero de canales: 3\n');
fprintf('Longitud de seal de entrada: %d muestras\n', L);
fprintf('Longitud de seal de salida: %d muestras\n', length(salidaMux));
fprintf('Frecuencia de salida: %d Hz\n', fs_salida);

for k = 1:3
    N_filtro = size(filtro_polifasico{k}, 2) * M;
    fprintf('Filtro %d - Orden total: %d, Subfiltros: %dx%d\n', ...
            k, N_filtro, M, size(filtro_polifasico{k}, 2));
end

%% Respuesta en frecuencia de los filtros
figure;
for k = 1:3
    % Reconstruir filtro original
    %h_original = reshape(filtro_polifasico{k}', 1, []);

    [H, w] = freqz(h_original{k}, 1, 1024, fs);
    
    subplot(3,1,k);
    plot(w, 20*log10(abs(H)));
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(sprintf('Respuesta en frecuencia - Filtro Kaiser %d', k));
    grid on;
    
    % Resaltar banda de paso
    hold on;
    ylims = ylim;
    line([bandas{k}(1) bandas{k}(1)], ylims, 'Color', 'r', 'LineStyle', '--');
    line([bandas{k}(2) bandas{k}(2)], ylims, 'Color', 'r', 'LineStyle', '--');
    ylim([-60, 5]);
    hold off;
end

%% COMPARACIN CON LA FUNCIN FILTER DE MATLAB
fprintf('\n=== COMPARACIN CON FUNCIN FILTER ===\n');
tic;


% Hacer upsampling de las seales (insertar 14 ceros entre cada muestra)
canal1_up = zeros(L*factor, 1);
canal2_up = zeros(L*factor, 1);
canal3_up = zeros(L*factor, 1);

multDesplazamintoFrecuencia_filter = 1;
for n = 1:L
    idx = (n-1)*factor + 1;
    
    % Aplicar desplazamiento en frecuencia
    muestra1 = canal1(n) * multDesplazamintoFrecuencia_filter;
    muestra2 = canal2(n);
    muestra3 = canal3(n) * multDesplazamintoFrecuencia_filter;
    multDesplazamintoFrecuencia_filter = multDesplazamintoFrecuencia_filter * (-1);
    
    % Upsampling: solo la primera muestra, las dems son cero
    canal1_up(idx) = muestra1;
    canal2_up(idx) = muestra2;
    canal3_up(idx) = muestra3;
end

% Filtrar usando la funcin filter()
canal1_filtrado = filter(h_original{1}, 1, canal1_up);
canal2_filtrado = filter(h_original{2}, 1, canal2_up);
canal3_filtrado = filter(h_original{3}, 1, canal3_up);

% Combinar canales
salidaMux_filter = 10*(canal1_filtrado + canal2_filtrado + canal3_filtrado);

tiempo_filter = toc;
fprintf('Tiempo de procesamiento con filter(): %.4f segundos\n', tiempo_filter);

%% ANLISIS COMPARATIVO
diferencia = salidaMux_filter - salidaMux;
error_rms = sqrt(mean(diferencia.^2));
error_max = max(abs(diferencia));
error_relativo = error_rms / sqrt(mean(salidaMux.^2));

fprintf('Error RMS entre mtodos: %.2e\n', error_rms);
fprintf('Error mximo absoluto: %.2e\n', error_max);
fprintf('Error relativo (%%): %.6f%%\n', error_relativo * 100);

%% COMPARACIN ESPECTRAL
N_comp = length(salidaMux_filter);
f_comp = linspace(0, fs_salida, N_comp);

X_filter = fft(salidaMux_filter);
X_polifasico = fft(salidaMux);

magX_filter_dB = 20*log10(abs(X_filter));
magX_polifasico_dB = 20*log10(abs(X_polifasico));

%% GRFICAS COMPARATIVAS
figure;
subplot(2,2,1);
plot(f_comp(1:N_comp/2)/1000, magX_filter_dB(1:N_comp/2), 'b-', 'LineWidth', 1.5);
hold on;
plot(f_comp(1:N_comp/2)/1000, magX_polifasico_dB(1:N_comp/2), 'r--', 'LineWidth', 1);
xlabel('Frecuencia (kHz)');
ylabel('Magnitud (dB)');
title('Comparacin: filter() vs Polifsico');
legend('filter()', 'Polifsico', 'Location', 'best');
grid on;

subplot(2,2,2);
t_comp = (0:min(2000, N_comp-1))/fs_salida;
plot(t_comp, salidaMux_filter(1:length(t_comp)), 'b-', 'LineWidth', 1.5);
hold on;
plot(t_comp, salidaMux(1:length(t_comp)), 'r--', 'LineWidth', 1);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Comparacin temporal');
legend('filter()', 'Polifsico');
grid on;

subplot(2,2,3);
plot(t_comp, diferencia(1:length(t_comp)), 'g-', 'LineWidth', 1);
xlabel('Tiempo (s)');
ylabel('Error');
title(sprintf('Error entre mtodos (RMS=%.2e)', error_rms));
grid on;

subplot(2,2,4);
diff_espectral = magX_filter_dB(1:N_comp/2) - magX_polifasico_dB(1:N_comp/2);
plot(f_comp(1:N_comp/2)/1000, diff_espectral, 'g-', 'LineWidth', 1);
xlabel('Frecuencia (kHz)');
ylabel('Diferencia (dB)');
title('Diferencia espectral');
grid on;

%% CONCLUSIONES
fprintf('\n=== RESULTADOS COMPARATIVOS ===\n');
if error_relativo < 1e-10
    fprintf('? Los mtodos son EQUIVALENTES (error < 1e-8%%)\n');
elseif error_relativo < 1e-6
    fprintf('? Los mtodos son PRCTICAMENTE IDNTICOS (error < 1e-4%%)\n');
else
    fprintf('? Los mtodos tienen diferencias (error = %.2e%%)\n', error_relativo*100);
end

save('comparacion_filter.mat', 'salidaMux_filter', 'diferencia', 'error_rms', 'tiempo_filter');

%% ANLISIS DE CANALES FILTRADOS POR SEPARADO
fprintf('\n=== ANLISIS DE CANALES FILTRADOS INDIVIDUALMENTE ===\n');

% Procesar cada canal por separado con la implementacin polifsica
canales_filtrados = cell(3,1);
nombres_canales = {'Canal 1 (12.3-15.4 kHz)', 'Canal 2 (16.3-19.4 kHz)', 'Canal 3 (20.3-23.4 kHz)'};
colores = {'b', 'r', 'g'};

for canal_idx = 1:3
    fprintf('Procesando %s...\n', nombres_canales{canal_idx});
    
    % Inicializar salida para este canal
    salida_canal = zeros(L*factor, 1);
    
    % Memoria para los subfiltros de este canal
    num_taps = size(filtro_polifasico{canal_idx}, 2);
    memoria_canal = zeros(factor, num_taps);
    
    multDesplazamintoFrecuencia = 1;
    
    % Seleccionar el canal de entrada correspondiente
    if canal_idx == 1
        canal_entrada = canal1;
    elseif canal_idx == 2
        canal_entrada = canal2;
    else
        canal_entrada = canal3;
    end
    
    % Procesar muestra por muestra
    for n = 1:L
        idx = (n-1)*factor + 1;
        
        % Aplicar desplazamiento solo a canales 1 y 3
        if canal_idx == 1 || canal_idx == 3
            muestra_actual = canal_entrada(n) * multDesplazamintoFrecuencia;
            multDesplazamintoFrecuencia = multDesplazamintoFrecuencia * (-1);
        else
            muestra_actual = canal_entrada(n);
        end
        
        % Inicializar salidas del filtrado
        salida_filtrada = zeros(factor, 1);
        
        % Procesar cada subfiltro polifsico
        H_poly = filtro_polifasico{canal_idx};
        
        for m = 1:factor
            % Desplazar memoria del subfiltro m
            memoria_canal(m, 2:end) = memoria_canal(m, 1:end-1);
            memoria_canal(m, 1) = muestra_actual;
            
            % Calcular salida del subfiltro m
            y_m = sum(H_poly(m, :) .* memoria_canal(m, :));
            salida_filtrada(m) = y_m;
        end
        
        % Guardar las M muestras de salida
        salida_canal(idx:idx+factor-1) = salida_filtrada;
    end
    
    % Guardar resultado
    canales_filtrados{canal_idx} = salida_canal;
    
    % Verificar que la suma coincida con la salida total
    if canal_idx == 1
        suma_verificacion = salida_canal;
    else
        suma_verificacion = suma_verificacion + salida_canal;
    end
end

% Verificar que la suma de canales individuales = salida total
error_verificacion = max(abs(suma_verificacion - salidaMux));
fprintf('Error de verificacin (suma canales vs total): %.2e\n', error_verificacion);

%% ANLISIS ESPECTRAL DE CADA CANAL
fprintf('\n--- Anlisis espectral individual ---\n');

N_sep = length(canales_filtrados{1});
f_sep = linspace(0, fs_salida, N_sep);
f_plot_sep = f_sep(1:N_sep/2);

% Calcular FFT de cada canal
espectros = cell(3,1);
for k = 1:3
    X_canal = fft(canales_filtrados{k});
    espectros{k} = 20*log10(abs(X_canal(1:N_sep/2)) + eps);
    
    % Calcular potencia en la banda objetivo
    if k == 1
        banda_obj = [12300, 15400];
    elseif k == 2
        banda_obj = [16300, 19400];
    else
        banda_obj = [20300, 23400];
    end
    
    idx_banda = find(f_plot_sep >= banda_obj(1) & f_plot_sep <= banda_obj(2));
    potencia_banda = sum(abs(X_canal(idx_banda)).^2);
    potencia_total = sum(abs(X_canal(1:N_sep/2)).^2);
    porcentaje_banda = 100 * potencia_banda / potencia_total;
    
    fprintf('Canal %d: %.1f%% de potencia en banda objetivo (%.1f-%.1f kHz)\n', ...
            k, porcentaje_banda, banda_obj(1)/1000, banda_obj(2)/1000);
end

%% GRFICAS DE CANALES INDIVIDUALES

% Grfica 1: Espectros individuales superpuestos
figure('Position', [100, 100, 1400, 900]);

subplot(3,2,1);
hold on;
for k = 1:3
    plot(f_plot_sep/1000, espectros{k}, colores{k}, 'LineWidth', 2, 'DisplayName', nombres_canales{k});
end

% Marcar bandas objetivo
bandas_obj = {[12.3 15.4], [16.3 19.4], [20.3 23.4]};
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

% Grfica 2: Zoom en cada banda
subplot(3,2,2);
hold on;
for k = 1:3
    % Zoom en la banda de cada canal
    f_zoom_start = bandas_obj{k}(1) - 2;
    f_zoom_end = bandas_obj{k}(2) + 2;
    idx_zoom = find(f_plot_sep/1000 >= f_zoom_start & f_plot_sep/1000 <= f_zoom_end);
    
    plot(f_plot_sep(idx_zoom)/1000, espectros{k}(idx_zoom), colores{k}, ...
         'LineWidth', 2, 'DisplayName', nombres_canales{k});
    
    % Marcar banda objetivo
    line([bandas_obj{k}(1) bandas_obj{k}(1)], ylim, 'Color', colores{k}, 'LineStyle', '--');
    line([bandas_obj{k}(2) bandas_obj{k}(2)], ylim, 'Color', colores{k}, 'LineStyle', '--');
end
xlabel('Frecuencia (kHz)');
ylabel('Magnitud (dB)');
title('Zoom en bandas objetivo');
legend('show', 'Location', 'best');
grid on;
hold off;

% Grficas 3-5: Cada canal por separado
for k = 1:3
    subplot(3,2,k+2);
    plot(f_plot_sep/1000, espectros{k}, colores{k}, 'LineWidth', 1.5);
    
    % Marcar banda objetivo
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
    
    % Aadir informacin de la banda
    text(bandas_obj{k}(1) + (bandas_obj{k}(2)-bandas_obj{k}(1))/2, ylims(2)-10, ...
         sprintf('%.1f-%.1f kHz', bandas_obj{k}(1), bandas_obj{k}(2)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', colores{k});
    hold off;
end

%% ANLISIS TEMPORAL DE CANALES
figure('Position', [200, 200, 1400, 600]);

t_temp = (0:min(5000, N_sep-1))/fs_salida;
for k = 1:3
    subplot(1,3,k);
    plot(t_temp, canales_filtrados{k}(1:length(t_temp)), colores{k}, 'LineWidth', 1);
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    title([nombres_canales{k} ' - Forma de onda']);
    grid on;
    
    % Estadsticas del canal
    rms_canal = sqrt(mean(canales_filtrados{k}.^2));
    max_canal = max(abs(canales_filtrados{k}));
    text(0.02, 0.95, sprintf('RMS: %.4f\nMx: %.4f', rms_canal, max_canal), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'white', 'EdgeColor', colores{k});
end

%% ESTADSTICAS DETALLADAS
fprintf('\n--- Estadsticas de canales filtrados ---\n');
for k = 1:3
    canal_data = canales_filtrados{k};
    
    fprintf('\n%s:\n', nombres_canales{k});
    fprintf('  RMS: %.6f\n', sqrt(mean(canal_data.^2)));
    fprintf('  Mximo: %.6f\n', max(abs(canal_data)));
    fprintf('  Potencia total: %.6f\n', sum(canal_data.^2)/length(canal_data));
    
    % Analizar contenido en diferentes bandas de frecuencia
    X_canal = fft(canal_data);
    
    % Banda objetivo
    banda_target = bandas_obj{k};
    idx_target = find(f_plot_sep >= banda_target(1)*1000 & f_plot_sep <= banda_target(2)*1000);
    potencia_target = sum(abs(X_canal(idx_target)).^2);
    
    % Fuera de banda (aliasing)
    idx_fuera = find(f_plot_sep < banda_target(1)*1000 | f_plot_sep > banda_target(2)*1000);
    potencia_fuera = sum(abs(X_canal(idx_fuera)).^2);
    
    relacion_dentro_fuera = 10*log10(potencia_target / potencia_fuera);
    fprintf('  Relacin dentro/fuera de banda: %.1f dB\n', relacion_dentro_fuera);
end

% Guardar canales individuales
save('canales_filtrados_individuales.mat', 'canales_filtrados', 'nombres_canales');
fprintf('\nCanales guardados en: canales_filtrados_individuales.mat\n');

