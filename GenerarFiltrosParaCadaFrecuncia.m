fs_list = [16000, 24000, 32000]; %Genero un vector con las fsamples (puedo hacerlo porque son conocidas)
% FILTROS DIGITALES PARA 3 TASAS DE MUESTREO
filtros = {};
fp = 3400;
Ap = 1;
As = 30;
ordenes = [];
for i = 1:3
    fs = fs_list(i);
    fprintf('\n--- Filtro digital para fs = %d Hz ---\n', fs);
    Wp = fp / (fs/2);  % Normalizaci�n de la Frecuencia de paso de 1dB a la frecuencia de sampleo de cada se�al
    [n, Wn] = buttord(Wp, 0.5, Ap, As);
    [b, a] = butter(n, Wn);
    filtros{i} = struct('fs', fs, 'b', b, 'a', a);
    ordenes(i) = n;
    % Graficar las respustas de cada apartado.
    figure;
    freqz(b, a, 1024, fs);
    title(sprintf('Filtro digital Butterworth - fs = %d Hz', fs));
    drawnow; 

    figure;
    zplane(b, a);
    title(sprintf('Polos y ceros - fs = %d Hz', fs));
    drawnow; 
    
    figure;
    impz(b, a, [], fs);
    title(sprintf('Respuesta al impulso - fs = %d Hz', fs));
    drawnow; 
end

save('filtros_guardados.mat', 'filtros', 'ordenes');
