%% Projet Télécom - Démo : Transmission de Texte (ASCII)
clear; clc; close all;

%% 1. Configuration
% --- Paramètres Utilisateur ---
message_texte = 'Bravo pour ce projet Telecom 2026 !'; % Ton message ici
SNR_dB = 1;  % Teste avec 4dB (bcp d'erreurs), 7dB (qqs erreurs), 12dB (parfait)
% ------------------------------

% Paramètres Système (QPSK)
M = 4; k = 2;
sps = 4; rolloff = 0.5; span = 6;
h = rcosdesign(rolloff, span, sps, 'sqrt'); % Filtre de mise en forme

%% 2. Émetteur (Transmitter)
fprintf('--- Émission ---\n');
fprintf('Message envoyé : "%s"\n', message_texte);

% A. Conversion Texte -> Bits
% 1. Convertir chaque char en nombre (uint8)
ascii_vals = uint8(message_texte);
% 2. Convertir en binaire (Matrice de char '0' et '1')
bin_matrix = dec2bin(ascii_vals, 8);
% 3. Convertir en vecteur colonne numérique
bits_temp = bin_matrix.' - '0'; % Transpose pour l'ordre des bits
dataIn = bits_temp(:);
numBits = length(dataIn);

% B. Modulation QPSK (Manuelle)
% Groupement par 2 bits
bits_reshaped = reshape(dataIn, k, [])';
sym_I = 2*bits_reshaped(:,1) - 1;
sym_Q = 2*bits_reshaped(:,2) - 1;
sym = (sym_I + 1j*sym_Q) / sqrt(2);

% C. Filtrage d'émission
txSignal = conv(upsample(sym, sps), h, 'same');

%% 3. Canal (Bruit AWGN)
SNR_lin = 10^(SNR_dB/10);
signalPower = mean(abs(txSignal).^2);
noisePower = signalPower / SNR_lin;

% Génération du bruit
noise = sqrt(noisePower/2) * (randn(size(txSignal)) + 1j*randn(size(txSignal)));
rxSignal = txSignal + noise;

%% 4. Récepteur (Receiver)
% A. Filtrage de réception + Sous-échantillonnage
rxFiltered = conv(rxSignal, h, 'same');
rxSym = rxFiltered(1:sps:end) / sum(h.^2); % Normalisation

% B. Démodulation QPSK (Décision rigide)
rx_bits_I = real(rxSym) > 0;
rx_bits_Q = imag(rxSym) > 0;

% Reconstitution du vecteur bits
dataOut_reshaped = [rx_bits_I rx_bits_Q];
dataOut = reshape(dataOut_reshaped', [], 1); % Vecteur colonne

%% 5. Reconstitution ASCII et Analyse
fprintf('\n--- Réception (SNR = %d dB) ---\n', SNR_dB);

% Calcul du nombre d'erreurs binaires
nb_erreurs_bits = sum(dataIn ~= dataOut);
BER = nb_erreurs_bits / numBits;

try
    % Conversion Bits -> Texte
    % 1. On reforme la matrice (8 lignes par caractère)
    nb_chars = length(dataOut) / 8;
    bits_recus_mat = reshape(dataOut, 8, nb_chars).';
    
    % 2. Conversion Matrice bits -> String binaire -> Décimal -> Char
    bin_str_rx = char(bits_recus_mat + '0');
    ascii_rx = bin2dec(bin_str_rx);
    message_recu = char(ascii_rx)'; % Transpose pour ligne
    
    % Affichage
    fprintf('Message reçu   : "%s"\n', message_recu);
    fprintf('Erreurs bits   : %d / %d (BER = %.4f)\n', nb_erreurs_bits, numBits, BER);
    
catch
    fprintf('Erreur critique : Impossible de reconstruire le texte.\n');
end

%% 6. Visualisation
figure('Name', 'Démo Transmission Texte', 'Position', [200, 200, 800, 400]);

% Constellation
subplot(1,2,1);
plot(real(rxSym), imag(rxSym), 'b.');
hold on;
plot(real(sym), imag(sym), 'r+', 'LineWidth', 2);
title(['Constellation Reçue (SNR = ' num2str(SNR_dB) ' dB)']);
xlabel('I'); ylabel('Q'); axis square; grid on;
legend('Reçu', 'Idéal');

% Comparaison Visuelle des bits (Zoom sur le début)
subplot(1,2,2);
L = min(100, length(dataIn)); % On affiche les 100 premiers bits
stem(1:L, dataIn(1:L), 'b', 'filled', 'DisplayName', 'Envoyé');
hold on;
stem(1:L, dataOut(1:L), 'r--', 'DisplayName', 'Reçu');
title('Comparaison bits (Zoom 100 premiers)');
ylim([-0.2 1.2]);
legend; grid on;