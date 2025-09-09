load('salidaMux.mat','salidaMux');

L = floor(length(salidaMux)/15);

salida1 = zeros(L,1);
salida2 = zeros(L,1);
salida3 = zeros(L,1);


for i = 1:L 
    %bloque = salidaMux(i*15 + 1 : (i+1)*15-1);  % vector 15x1

    idx_inicio = (i-1)*15 + 1;
    idx_fin    = idx_inicio + 15 - 1;
    bloque = salidaMux(idx_inicio:idx_fin);  % vector 15x1
    salida1(i) = dot(bloque, filtro{1}); 
    salida2(i) = dot(bloque, filtro{2}); 
    salida3(i) = dot(bloque, filtro{3}); 
    
end
load('salidaMux.mat','salidaMux');

% Número de bloques de 15 muestras
L = floor(length(salidaMux)/15);

% Inicializar vectores de salida
salida1 = zeros(L,1);
salida2 = zeros(L,1);
salida3 = zeros(L,1);
multDesplazamintoFrecuncia = 1;
for i = 1:length(filtro{1})
    filtro{1}(i) = filtro{1}(i) * multDesplazamintoFrecuncia;
    filtro{3}(i) = filtro{3}(i) * multDesplazamintoFrecuncia;
    multDesplazamintoFrecuncia = multDesplazamintoFrecuncia * (-1);
end

% Demultiplexado usando FIR
for i = 1:L 
    idx_inicio = (i-1)*15 + 1;
    idx_fin    = idx_inicio + 15 - 1;
    bloque = salidaMux(idx_inicio:idx_fin);  % vector 15x1

    salida1(i) = 10*dot(bloque, filtro{1}); 
    salida2(i) = 10*dot(bloque, filtro{2}); 
    salida3(i) = 10*dot(bloque, filtro{3}); 
end

% Guardar las salidas demultiplexadas en .mat
save('salidas_demultiplexadas.mat','salida1','salida2','salida3');

% Frecuencia de muestreo
fs_salida = 8000;  % Hz

% Lista de salidas y nombres de archivos
salidas = {salida1, salida2, salida3};
titulos = {'Canal 1','Canal 2','Canal 3'};
archivos = {'canal1_salida.wav','canal2_salida.wav','canal3_salida.wav'};

for k = 1:3
    x = salidas{k};

    % Guardar archivo WAV
    audiowrite(archivos{k}, x, fs_salida);
    fprintf('Archivo guardado: %s\n', archivos{k});

    % Reproducir audio
    fprintf('Reproduciendo %s...\n', titulos{k});
    sound(x, fs_salida);
    pause(length(x)/fs_salida + 0.5);  % Esperar hasta que termine la reproducción

    % FFT y espectros
    N = length(x);
    X = fft(x);
    f = linspace(0, fs_salida, N);
    magX = 20*log10(abs(X));

    % Graficar espectro completo
    figure;
    plot(f, magX);
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(['Espectro completo - ' titulos{k}]);
    grid on;

    % Graficar espectro hasta Nyquist
    figure;
    plot(f(1:N/2), magX(1:N/2));
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(['Espectro hasta Nyquist - ' titulos{k}]);
    grid on;

    % Graficar forma de onda
    t = (0:N-1)/fs_salida;
    figure;
    plot(t, x);
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    title(['Forma de onda - ' titulos{k}]);
    grid on;
end