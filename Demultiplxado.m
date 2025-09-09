%% DEMULTIPLEXOR FDM POLIFÁSICO EFICIENTE - VERSIÓN CORREGIDA
clear; clc;

% Cargar señal multiplexada y filtros
load('salidaMux_polifasico.mat','salidaMux');
load('filtro_polifasico.mat', 'filtro_polifasico');

fprintf('=== DEMULTIPLEXOR FDM POLIFÁSICO EFICIENTE ===\n');

%% PARÁMETROS DEL SISTEMA
factor_up = 15;           
factor_down = 15;         
fs_entrada = 120000;      
fs_salida = 8000;         

N_entrada = length(salidaMux);
L_salida = floor(N_entrada / factor_down);

fprintf('Señal de entrada: %d muestras a %d Hz\n', N_entrada, fs_entrada);
fprintf('Señal de salida: %d muestras a %d Hz por canal\n', L_salida, fs_salida);

%% PREPARAR FILTROS DEMULTIPLEXORES
filtros_demux = cell(3,1);
for k = 1:3
    H_poly = filtro_polifasico{k};
    filtros_demux{k} = reshape(H_poly', 1, []);
    fprintf('Filtro demux %d: Orden %d\n', k, length(filtros_demux{k})-1);
end

% Aplicar desplazamiento en frecuencia a filtros 1 y 3
fprintf('Aplicando desplazamiento en frecuencia...\n');
for k = [1, 3]
    multDesplazamiento = 1;
    for i = 1:length(filtros_demux{k})
        filtros_demux{k}(i) = filtros_demux{k}(i) * multDesplazamiento;
        multDesplazamiento = multDesplazamiento * (-1);
    end
end

%% IMPLEMENTACIÓN EFICIENTE DEL DEMULTIPLEXOR
fprintf('Iniciando demultiplexado eficiente...\n');
tic;

% Inicializar salidas
salida1 = zeros(L_salida, 1);
salida2 = zeros(L_salida, 1);
salida3 = zeros(L_salida, 1);

% Memorias independientes para cada filtro demultiplexor
memoria1 = zeros(length(filtros_demux{1}), 1);
memoria2 = zeros(length(filtros_demux{2}), 1);
memoria3 = zeros(length(filtros_demux{3}), 1);

% Índices independientes para las memorias circulares
idx_memoria1 = 1;
idx_memoria2 = 1;
idx_memoria3 = 1;

fprintf('Procesando %d muestras de entrada...\n', N_entrada);

for n = 1:N_entrada
    if mod(n, 50000) == 0
        fprintf('Procesada muestra %d de %d (%.1f%%)\n', n, N_entrada, 100*n/N_entrada);
    end
    
    % Actualizar memorias circulares con la nueva muestra (cada canal independiente)
    memoria1(idx_memoria1) = salidaMux(n);
    memoria2(idx_memoria2) = salidaMux(n);
    memoria3(idx_memoria3) = salidaMux(n);
    
    % Solo calcular salida cada "factor_down" muestras (downsampling eficiente)
    if mod(n, factor_down) == 0
        idx_salida = n / factor_down;
        
        % Canal 1
        suma1 = 0;
        for tap = 1:length(filtros_demux{1})
            mem_idx = mod(idx_memoria1 - tap, length(filtros_demux{1})) + 1;
            suma1 = suma1 + filtros_demux{1}(tap) * memoria1(mem_idx);
        end
        salida1(idx_salida) = suma1;
        
        % Canal 2  
        suma2 = 0;
        for tap = 1:length(filtros_demux{2})
            mem_idx = mod(idx_memoria2 - tap, length(filtros_demux{2})) + 1;
            suma2 = suma2 + filtros_demux{2}(tap) * memoria2(mem_idx);
        end
        salida2(idx_salida) = suma2;
        
        % Canal 3
        suma3 = 0;
        for tap = 1:length(filtros_demux{3})
            mem_idx = mod(idx_memoria3 - tap, length(filtros_demux{3})) + 1;
            suma3 = suma3 + filtros_demux{3}(tap) * memoria3(mem_idx);
        end
        salida3(idx_salida) = suma3;
    end
    
    % Actualizar índices de memoria circular independientes
    idx_memoria1 = mod(idx_memoria1, length(filtros_demux{1})) + 1;
    idx_memoria2 = mod(idx_memoria2, length(filtros_demux{2})) + 1;
    idx_memoria3 = mod(idx_memoria3, length(filtros_demux{3})) + 1;
end

tiempo_demux = toc;
fprintf('Demultiplexado completado en %.4f segundos\n', tiempo_demux);

% Factor de ganancia para compensar el procesamiento
factor_ganancia = 10;
salida1 = salida1 * factor_ganancia;
salida2 = salida2 * factor_ganancia;
salida3 = salida3 * factor_ganancia;

fprintf('Aplicado factor de ganancia: %d\n', factor_ganancia);

%% GUARDAR RESULTADOS
save('salidas_demultiplexadas_eficientes.mat', 'salida1', 'salida2', 'salida3');

%% ANÁLISIS Y REPRODUCCIÓN
salidas = {salida1, salida2, salida3};
titulos = {'Canal 1 Demultiplexado', 'Canal 2 Demultiplexado', 'Canal 3 Demultiplexado'};
archivos_wav = {'canal1_demux.wav', 'canal2_demux.wav', 'canal3_demux.wav'};
colores = {'b', 'r', 'g'};

for k = 1:3
    fprintf('\n--- %s ---\n', titulos{k});
    x = salidas{k};
    
    % Estadísticas
    rms_val = sqrt(mean(x.^2));
    max_val = max(abs(x));
    fprintf('RMS: %.6f, Máximo: %.6f\n', rms_val, max_val);
    
    % Guardar WAV
    x_norm = x / max(abs(x)) * 0.95;
    audiowrite(archivos_wav{k}, x_norm, fs_salida);
    fprintf('Audio guardado: %s\n', archivos_wav{k});
    
    % Análisis espectral
    N = length(x);
    X = fft(x);
    f = linspace(0, fs_salida, N);
    [~, idx_max] = max(abs(X(2:N/2)));
    freq_dominante = f(idx_max + 1);
    fprintf('Frecuencia dominante: %.1f Hz\n', freq_dominante);
end

%% GRÁFICAS
figure('Position', [100, 100, 1400, 900]);

for k = 1:3
    x = salidas{k};
    N = length(x);
    X = fft(x);
    f = linspace(0, fs_salida, N);
    magX_dB = 20*log10(abs(X) + eps);
    
    % Espectro completo
    subplot(3,2,k*2-1);
    plot(f, magX_dB, colores{k}, 'LineWidth', 1.5);
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title([titulos{k} ' - Espectro completo']);
    grid on;
    
    % Espectro hasta Nyquist
    subplot(3,2,k*2);
    plot(f(1:N/2), magX_dB(1:N/2), colores{k}, 'LineWidth', 1.5);
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title([titulos{k} ' - Hasta Nyquist']);
    grid on;
end

% Formas de onda
figure('Position', [200, 200, 1400, 600]);
for k = 1:3
    x = salidas{k};
    t = (0:length(x)-1)/fs_salida;
    
    subplot(1,3,k);
    plot(t, x, colores{k}, 'LineWidth', 1);
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    title([titulos{k} ' - Forma de onda']);
    grid on;
    
    if length(t) > 2*fs_salida
        xlim([0, 2]);
    end
end

% Comparación espectral
figure('Position', [300, 300, 1200, 600]);
subplot(1,2,1);
hold on;
for k = 1:3
    x = salidas{k};
    N = length(x);
    X = fft(x);
    f = linspace(0, fs_salida, N);
    magX_dB = 20*log10(abs(X) + eps);
    plot(f(1:N/2), magX_dB(1:N/2), colores{k}, 'LineWidth', 2, 'DisplayName', titulos{k});
end
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Comparación espectral');
legend('show');
grid on;
hold off;

subplot(1,2,2);
hold on;
for k = 1:3
    x = salidas{k};
    t = (0:min(1000, length(x)-1))/fs_salida;
    plot(t, x(1:length(t)), colores{k}, 'LineWidth', 1.5, 'DisplayName', titulos{k});
end
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Comparación temporal');
legend('show');
grid on;
hold off;

fprintf('\n=== DEMULTIPLEXADO COMPLETADO ===\n');