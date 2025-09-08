
load('filtros_guardados.mat', 'filtros');



%% FILTRADO Y ESPECTRO DE SALIDA
fprintf('\n--- Filtrado digital de las se�ales ---\n');
for i = 1:3
    fdata = filtros{i};
for k = 1:3
    [x, fs] = audioread(archivos{i}{k});
    y = filter(fdata.b, fdata.a, x); %Esta notacion toma los elementos del struct despues del punto
    % Espectro
    Y = abs(fft(y));
    N = length(Y);
    f = linspace(0, fs, N);
    figure(i);
    subplot(3,1,k);
    plot(f(1:N/2), 20*log10(Y(1:N/2)), 'Color', colores(k,:));
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(['Filtrada: ' archivos{i}{k}]);
    grid on;
end
end