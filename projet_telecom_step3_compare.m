%% Projet Télécom - Étape 3 : Comparaison QPSK vs 16-QAM
clear; clc; close all;

%% 1. Paramètres
numBits = 200000;       % Grand nombre de bits pour lisser les courbes
SNR_dB_range = 0:1:18;  % On va plus loin en SNR car le 16-QAM a besoin de plus de puissance

% Stockage des résultats
ber_qpsk = zeros(size(SNR_dB_range));
ber_16qam = zeros(size(SNR_dB_range));

%% --- SIMULATION QPSK (Rappel rapide) ---
fprintf('Simulation QPSK en cours...\n');
M = 4; k = 2;
nBits_Q = ceil(numBits/k)*k;
dataIn_Q = randi([0 1], nBits_Q, 1);

% Mod QPSK (Mapping simple)
dataReshaped = reshape(dataIn_Q, k, length(dataIn_Q)/k)';
sym_I = 2*dataReshaped(:,1) - 1;
sym_Q = 2*dataReshaped(:,2) - 1;
txSig_Q = (sym_I + 1j*sym_Q) / sqrt(2); % Normalisation P=1

%% --- SIMULATION 16-QAM (Nouvelle partie) ---
fprintf('Simulation 16-QAM en cours...\n');
M16 = 16; k16 = 4;
nBits_16 = ceil(numBits/k16)*k16;
dataIn_16 = randi([0 1], nBits_16, 1);

% Mod 16-QAM Manuelle
% On groupe par 4 bits : 2 pour I, 2 pour Q
dataReshaped16 = reshape(dataIn_16, k16, length(dataIn_16)/k16)';
bits_I_16 = dataReshaped16(:, 1:2); % 2 premiers bits
bits_Q_16 = dataReshaped16(:, 3:4); % 2 derniers bits

% Fonction locale pour convertir 2 bits en niveau PAM4 (-3, -1, 1, 3)
% Bit à bit -> Decimal: 00->0, 01->1, 10->2, 11->3
dec_I = bits_I_16(:,1)*2 + bits_I_16(:,2); 
dec_Q = bits_Q_16(:,1)*2 + bits_Q_16(:,2);

% Mapping vers niveaux: 0->-3, 1->-1, 2->1, 3->3
% Formule mathématique : niveau = 2*decimal - 3
lvl_I = 2*dec_I - 3;
lvl_Q = 2*dec_Q - 3;

txSig_16 = lvl_I + 1j*lvl_Q;

% Normalisation de puissance 16-QAM
% Puissance moyenne théorique = (3^2 + 1^2 + (-1)^2 + (-3)^2) / 4 = 5 par dimension
% P_total = 10. Donc on divise par sqrt(10).
txSig_16 = txSig_16 / sqrt(10);

%% 3. Boucle SNR (Commune)
for i = 1:length(SNR_dB_range)
    SNR_dB = SNR_dB_range(i);
    SNR_lin = 10^(SNR_dB/10);
    noisePower = 1 / SNR_lin;
    noiseScale = sqrt(noisePower/2);
    
    % --- CANAL QPSK ---
    noise_Q = noiseScale * (randn(size(txSig_Q)) + 1j*randn(size(txSig_Q)));
    rxSig_Q = txSig_Q + noise_Q;
    
    % Démodulation QPSK
    rx_demod_Q = rxSig_Q * sqrt(2);
    rx_bits_I = real(rx_demod_Q) > 0;
    rx_bits_Q = imag(rx_demod_Q) > 0;
    dataOut_Q = reshape([rx_bits_I rx_bits_Q]', nBits_Q, 1);
    ber_qpsk(i) = sum(dataIn_Q ~= dataOut_Q) / nBits_Q;
    
    % --- CANAL 16-QAM ---
    noise_16 = noiseScale * (randn(size(txSig_16)) + 1j*randn(size(txSig_16)));
    rxSig_16 = txSig_16 + noise_16;
    
    % Démodulation 16-QAM
    rx_demod_16 = rxSig_16 * sqrt(10); % Dénormalisation
    
    % Décision seuils pour 16-QAM (-3, -1, 1, 3)
    % Seuils optimaux : -2, 0, 2
    
    % Décodage I
    rx_I_raw = real(rx_demod_16);
    dec_I_RX = zeros(size(rx_I_raw));
    dec_I_RX(rx_I_raw < -2) = 0;              % -> -3 (bits 00)
    dec_I_RX(rx_I_raw >= -2 & rx_I_raw < 0) = 1; % -> -1 (bits 01)
    dec_I_RX(rx_I_raw >= 0 & rx_I_raw < 2) = 2;  % -> +1 (bits 10)
    dec_I_RX(rx_I_raw >= 2) = 3;              % -> +3 (bits 11)
    
    % Décodage Q (identique)
    rx_Q_raw = imag(rx_demod_16);
    dec_Q_RX = zeros(size(rx_Q_raw));
    dec_Q_RX(rx_Q_raw < -2) = 0;
    dec_Q_RX(rx_Q_raw >= -2 & rx_Q_raw < 0) = 1;
    dec_Q_RX(rx_Q_raw >= 0 & rx_Q_raw < 2) = 2;
    dec_Q_RX(rx_Q_raw >= 2) = 3;
    
    % Reconstitution bits (Decimal -> Binaire manuel)
    % Bit 1 (MSB) = floor(dec / 2)
    % Bit 2 (LSB) = mod(dec, 2)
    b1_I = floor(dec_I_RX/2); b2_I = mod(dec_I_RX, 2);
    b1_Q = floor(dec_Q_RX/2); b2_Q = mod(dec_Q_RX, 2);
    
    dataOut_16 = reshape([b1_I b2_I b1_Q b2_Q]', nBits_16, 1);
    ber_16qam(i) = sum(dataIn_16 ~= dataOut_16) / nBits_16;
end

%% 4. Visualisation Comparée
figure;
semilogy(SNR_dB_range, ber_qpsk, 'b-o', 'LineWidth', 2, 'DisplayName', 'QPSK (Simulé)');
hold on;
semilogy(SNR_dB_range, ber_16qam, 'r-s', 'LineWidth', 2, 'DisplayName', '16-QAM (Simulé)');
grid on;
xlabel('SNR (dB)');
ylabel('BER (Bit Error Rate)');
title('Comparaison BER : QPSK vs 16-QAM');
legend;
ylim([1e-5 1]);

% Affichage d'une constellation 16-QAM bruitée (pour le plaisir)
figure;
plot(real(rxSig_16(1:2000)), imag(rxSig_16(1:2000)), 'r.');
hold on;
plot(real(txSig_16(1:100)), imag(txSig_16(1:100)), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 5);
title('Constellation 16-QAM (SNR élevé)');
grid on; axis square;