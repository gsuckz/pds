%% FILTRO ANTIALIASING ANAL�GICO (para fs = 8 kHz)
%% Dise�o del filtro analogico necesario (No implementado, solo se grafica su funci�n de Transferencia)
fp = 3400; Ap = 1;
fs_alias = 4000; As = 40;

fprintf('\n--- Filtro antialiasing analogico ---\n')
[na, Wn_a] = buttord(2*pi*fp, 2*pi*fs_alias, Ap, As, 's');
[ba, aa] = butter(na, Wn_a, 's');
Ha = tf(ba, aa);

figure;
bode(Ha);
title('Filtro antialiasing analogico');
grid on;