%% Cargar archivos procesados a 8 kHz
load('archivos_procesados.mat', 'archivos_procesados');
colores = lines(3);

fs_in = 8000;      % fs original
fs_out = 120000;   % fs deseada
factor = fs_out / fs_in;  % debe ser entero (120k / 8k = 15)

if mod(factor,1) ~= 0
    error('El factor de interpolación no es entero.');
end

fprintf('\n--- Generando señales a 120 kHz ---\n');
for k = 1:3
    % Leer archivo procesado (8 kHz)
    [x8k, fs] = audioread(archivos_procesados{k});
    if fs ~= fs_in
        error('El archivo no tiene fs = 8 kHz como se esperaba.');
    end
    
    % Interpolación (upsample) a 120 kHz
    x120k = resample(x8k, fs_out, fs_in);
    
    % Guardar archivo a 120 kHz
    [filepath, name, ext] = fileparts(archivos_procesados{k});
    nombre_salida = fullfile(filepath, [name '_120kHz' ext]);
    audiowrite(nombre_salida, x120k, fs_out);
    fprintf('Archivo guardado: %s\n', nombre_salida);

    % Calcular espectro
    Y = abs(fft(x120k));
    N = length(Y);
    f = linspace(0, fs_out, N);

    figure(100);
    subplot(3,1,k);
    plot(f, 20*log10(Y + eps), 'Color', colores(k,:));
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(sprintf('Canal %d procesado a 120 kHz', k));
    grid on;
    xlim([0 fs_out/2]);  % hasta Nyquist
end
