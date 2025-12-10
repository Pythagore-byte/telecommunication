%% Projet Télécom - Étape 1 : Chaîne de transmission de base (Sans Toolbox)
clear; clc; close all;

%% 1. Paramètres de simulation
M = 4;                  % QPSK (4 points)
k = log2(M);            % 2 bits par symbole
numBits = 10000;        % Nombre de bits
numBits = ceil(numBits/k)*k; % Ajustement taille

%% 2. Émetteur (Transmitter TX)
% Génération aléatoire de bits (0 ou 1)
dataIn = randi([0 1], numBits, 1); 

% --- MODULATION MANUELLE (QPSK) ---
% On regroupe les bits par paquets de 2 (I et Q)
dataReshaped = reshape(dataIn, k, length(dataIn)/k)';
bits_I = dataReshaped(:, 1); % Premier bit de chaque paire
bits_Q = dataReshaped(:, 2); % Second bit de chaque paire

% Mappage : 0 -> -1, 1 -> +1
sym_I = 2*bits_I - 1; 
sym_Q = 2*bits_Q - 1;  

% Création du signal complexe
dataMod = sym_I + 1j*sym_Q;
% Normalisation pour avoir une puissance moyenne de 1
% Pour QPSK non normalisé, la puissance est I^2 + Q^2 = 1^2 + 1^2 = 2.
% Il faut donc diviser par sqrt(2).
dataMod = dataMod / sqrt(2);

%% 3. Canal (Idéal)
receivedSignal = dataMod;

%% 4. Récepteur (Receiver RX)
% --- DÉMODULATION MANUELLE ---
% On enlève la normalisation pour retrouver les niveaux +/- 1
rx_demod_scaled = receivedSignal * sqrt(2);

% Décision sur la partie réelle (I)
rx_bits_I = real(rx_demod_scaled) > 0; % Si > 0 alors bit est 1, sinon 0

% Décision sur la partie imaginaire (Q)
rx_bits_Q = imag(rx_demod_scaled) > 0; % Si > 0 alors bit est 1, sinon 0

% Reconstitution du flux binaire
dataOutMatrix = [rx_bits_I, rx_bits_Q];
dataOut = reshape(dataOutMatrix', numBits, 1);

%% 5. Analyse et Visualisation

%% 5. Analyse et Visualisation

% Pour cette première étape : 1 échantillon par symbole
sps = 1;
txSignal = dataMod;

% --- A. Domaine Temporel ---
figure('Name', 'Analyses Temporelle et Fréquentielle');
subplot(2,1,1);

Nplot = min(200, length(txSignal));       % on trace max 200 échantillons
t_axis = (0:Nplot-1)/sps;                % temps en "durées symbole"

plot(t_axis, real(txSignal(1:Nplot)), 'b-', 'LineWidth', 1.5); hold on;
plot(t_axis, imag(txSignal(1:Nplot)), 'r--', 'LineWidth', 1.5);
grid on;
title('Signal Transmis (Domaine Temporel - Zoom)');
xlabel('Temps (en durées symbole)');
ylabel('Amplitude');
legend('Voie I (In-Phase)', 'Voie Q (Quadrature)');

% --- B. Domaine Fréquentiel (Densité Spectrale de Puissance) ---
subplot(2,1,2);
L = length(txSignal);
nfft = 2^nextpow2(L);
Y = fft(txSignal, nfft);
f = (-nfft/2:nfft/2-1)/(nfft/sps);   % fréquence normalisée à la fréquence symbole

power_spectrum = fftshift(abs(Y).^2/L);
power_spectrum_dB = 10*log10(power_spectrum + eps);  % +eps pour éviter log(0)

plot(f, power_spectrum_dB, 'k', 'LineWidth', 1.5);
grid on;
title('Densité Spectrale de Puissance (PSD)');
xlabel('Fréquence Normalisée (f / Rs)');
ylabel('Puissance (dB)');
xlim([-2 2]);
ylim([max(power_spectrum_dB)-60, max(power_spectrum_dB)+5]);

% a) Calcul du BER (Bit Error Rate) manuellement
errors = sum(dataIn ~= dataOut);
bit_error_rate = errors / numBits;

fprintf('--- Résultats de la simulation ---\n');
fprintf('Modulation : QPSK (Manuelle)\n');
fprintf('Nombre d''erreurs : %d\n', errors);
fprintf('Taux d''erreur binaire (BER) : %.4e\n', bit_error_rate);

% b) Visualisation
figure;
plot(real(receivedSignal), imag(receivedSignal), 'b.');
title('Constellation reçue (QPSK sans bruit)');
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
axis([-1.5 1.5 -1.5 1.5]);
grid on;
axis square;