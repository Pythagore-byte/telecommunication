%% Projet Télécom - Étape 1-Bis : Chaîne Complète avec Filtrage et Spectre
clear; clc; close all;

%% 1. Paramètres du système
% Je définis les paramètres de modulation et de filtrage
M = 4;                  % Je choisis la QPSK (4 états)
k = log2(M);            % Je calcule le nombre de bits par symbole (2 bits)
numBits = 10000;        % Je définis le nombre de bits à transmettre
numBits = ceil(numBits/k)*k; % Je m'assure que le nombre de bits est un multiple de k

% Paramètres du filtre de mise en forme (Pulse Shaping)
sps = 8;                % Samples Per Symbol: Je prends 8 échantillons par symbole pour avoir une bonne résolution
rolloff = 0.5;          % Facteur de retombée (beta): Je choisis 0.5 pour un compromis bande/temps
span = 6;               % Je définis la longueur d

%% 2. Émetteur (Transmitter TX)du filtre en nombre de symboles
% --- A. Génération des données ---
% Je génère une séquence binaire aléatoire
dataIn = randi([0 1], numBits, 1); %//ici numbits ligne et une colonne

% --- B. Modulation en Bande de Base (QPSK Manuelle) ---
% Je transforme les bits en symboles complexes
dataReshaped = reshape(dataIn, k, length(dataIn)/k)';
bits_I = dataReshaped(:, 1);
bits_Q = dataReshaped(:, 2);
sym_I = 2*bits_I - 1; 
sym_Q = 2*bits_Q - 1;
symbols = (sym_I + 1j*sym_Q) / sqrt(2); % Je normalise la puissance

% --- C. Sur-échantillonnage (Upsampling) ---
% J'insère des zéros entre les symboles pour préparer le filtrage numérique
symbols_upsampled = zeros(length(symbols)*sps, 1);
symbols_upsampled(1:sps:end) = symbols;

% --- D. Filtrage de mise en forme (RRC - Root Raised Cosine) ---
% Je génère la réponse impulsionnelle du filtre RRC
h = rcosdesign(rolloff, span, sps, 'sqrt'); 

% Je réalise la convolution pour obtenir le signal temporel "continu"
txSignal = conv(symbols_upsampled, h, 'same');

%% 3. Analyse Spectrale et Temporelle (Visualisation)

% --- A. Domaine Temporel ---
figure('Name', 'Analyses Temporelle et Fréquentielle');
subplot(2,1,1);
t_axis = (0:200-1)/sps; % Axe temps normalisé en symboles
plot(t_axis, real(txSignal(1:200)), 'b-', 'LineWidth', 1.5); hold on;
plot(t_axis, imag(txSignal(1:200)), 'r--', 'LineWidth', 1.5);
grid on;
title('Signal Transmis (Domaine Temporel - Zoom)');
xlabel('Temps (en durée symbole)');
ylabel('Amplitude');
legend('Voie I (In-Phase)', 'Voie Q (Quadrature)');
% Commentaire: On voit bien ici que le signal n'est plus "carré" mais lisse.

% --- B. Domaine Fréquentiel (Densité Spectrale de Puissance) ---
subplot(2,1,2);
L = length(txSignal);
nfft = 2^nextpow2(L); % Je prends une puissance de 2 pour optimiser la FFT
Y = fft(txSignal, nfft);
f = (-nfft/2:nfft/2-1)/(nfft/sps); % Axe fréquence normalisé
power_spectrum = fftshift(abs(Y).^2/L);
power_spectrum_dB = 10*log10(power_spectrum);

plot(f, power_spectrum_dB, 'k', 'LineWidth', 1.5);
grid on;
title('Densité Spectrale de Puissance (PSD)');
xlabel('Fréquence Normalisée (f / Rs)');
ylabel('Puissance (dB)');
xlim([-2 2]); ylim([max(power_spectrum_dB)-60 max(power_spectrum_dB)+5]);
% Commentaire: Je vérifie que l'énergie est bien confinée dans la bande passante.

%% 4. Récepteur (RX) - Filtrage Adapté et Décision

% --- A. Filtrage Adapté (Matched Filter) ---
% J'applique le même filtre en réception pour maximiser le SNR
rxSignal_filtered = conv(txSignal, h, 'same');

% --- B. Sous-échantillonnage (Downsampling) ---
% Je récupère un échantillon tous les 'sps' instants (l'instant optimal)
% Je dois tenir compte du délai introduit par le filtre (span*sps / 2)
delay = span*sps/2; 
% Note: 'conv' avec 'same' centre le signal, donc on échantillonne au milieu
rxSymbols = rxSignal_filtered(1:sps:end); 

% Normalisation du gain introduit par le filtre
rxSymbols = rxSymbols / (sum(h.^2)); 
% Petite correction d'amplitude due à la chaîne de filtrage
rxSymbols = rxSymbols / sqrt(10); % Ajustement empirique selon l'énergie du filtre

% --- C. Démodulation et BER ---
% Décision sur le signe (QPSK)
rx_bits_I = real(rxSymbols) > 0;
rx_bits_Q = imag(rxSymbols) > 0;
dataOutMatrix = [rx_bits_I, rx_bits_Q];
dataOut = reshape(dataOutMatrix', numBits, 1);

% Calcul des erreurs
errors = sum(dataIn ~= dataOut);
BER = errors / numBits;

fprintf('--- Validation Étape 1-Bis ---\n');
fprintf('BER après filtrage Tx/Rx : %.4f\n', BER);
fprintf('Si le BER est ~0, la chaîne de filtrage est correcte.\n');

% Visualisation de la constellation après filtrage adapté
figure;
plot(real(rxSymbols), imag(rxSymbols), 'b.');
title('Constellation après Filtrage Adapté (Sans Bruit)');
grid on; axis square;