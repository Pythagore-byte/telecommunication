%% oui
// %% Projet Télécom - Étape 2 : Ajout de Bruit (AWGN) et Courbes BER
// clear; clc; close all;

// %% 1. Paramètres
// numBits = 100000;       % On augmente le nombre de bits pour avoir des stats fiables
// M = 4;                  % QPSK
// k = log2(M);
// numBits = ceil(numBits/k)*k;

// SNR_dB_range = 0:1:12;  % Plage de SNR en dB à tester
// ber_simule = zeros(size(SNR_dB_range)); % Stockage des résultats

// %% 2. Émetteur (TX) - Fixe
// dataIn = randi([0 1], numBits, 1);
// dataReshaped = reshape(dataIn, k, length(dataIn)/k)';
// bits_I = dataReshaped(:, 1);
// bits_Q = dataReshaped(:, 2);

// % Modulation QPSK (Normalisée à 1 Watt)
// sym_I = (2*bits_I - 1);
// sym_Q = (2*bits_Q - 1);
// txSig = (sym_I + 1j*sym_Q) / sqrt(2);

// %% 3. Boucle de Simulation sur le SNR
// fprintf('Simulation en cours...\n');

// for i = 1:length(SNR_dB_range)
//     SNR_dB = SNR_dB_range(i);
    
//     % --- AJOUT DU BRUIT (MANUEL) ---
//     % 1. Convertir SNR dB en échelle linéaire
//     % SNR = P_signal / P_noise
//     % Comme P_signal = 1, alors P_noise = 1 / SNR_lineaire
//     SNR_lin = 10^(SNR_dB/10);
//     noisePower = 1 / SNR_lin;
    
//     % 2. Générer le bruit complexe
//     % randn génère un bruit de variance 1. On multiplie par sqrt(Puissance)
//     % On divise la puissance du bruit par 2 car elle se répartit sur I et Q
//     noiseScale = sqrt(noisePower/2);
//     noise = noiseScale * (randn(size(txSig)) + 1j*randn(size(txSig)));
    
//     % 3. Signal reçu
//     rxSig = txSig + noise;
    
//     % --- RÉCEPTEUR (RX) ---
//     % Démodulation (décision rigide)
//     rx_demod_scaled = rxSig * sqrt(2); % Dénormalisation
//     rx_bits_I = real(rx_demod_scaled) > 0;
//     rx_bits_Q = imag(rx_demod_scaled) > 0;
    
//     % Reconstitution
//     dataOutMatrix = [rx_bits_I, rx_bits_Q];
//     dataOut = reshape(dataOutMatrix', numBits, 1);
    
//     % Calcul des erreurs
//     errors = sum(dataIn ~= dataOut);
//     ber_simule(i) = errors / numBits;
// end

// %% 4. Comparaison Théorique (QPSK)
// % Formule théorique approximative pour QPSK avec mappage de Gray
// % BER = 0.5 * erfc(sqrt(Eb/N0))
// % Pour QPSK: SNR = 2 * Eb/N0  => Eb/N0 = SNR_lin / 2
// SNR_lin_range = 10.^(SNR_dB_range/10);
// ber_theorique = 0.5 * erfc(sqrt(SNR_lin_range / 2)); 

// %% 5. Visualisation
// figure;
// semilogy(SNR_dB_range, ber_simule, 'bo-', 'LineWidth', 2, 'DisplayName', 'Simulé');
// hold on;
// semilogy(SNR_dB_range, ber_theorique, 'r--', 'LineWidth', 2, 'DisplayName', 'Théorique');
// grid on;
// xlabel('SNR (dB)');
// ylabel('BER (Bit Error Rate)');
// title('Performance BER de la QPSK sur canal AWGN');
// legend;
// axis([0 12 1e-5 1]);

// % Visualisation d'une constellation bruitée (pour un SNR moyen, ex: 10dB)
// figure;
// SNR_demo = 10; 
// SNR_lin_demo = 10^(SNR_demo/10);
// noiseScale_demo = sqrt((1/SNR_lin_demo)/2);
// rxSig_demo = txSig + noiseScale_demo * (randn(size(txSig)) + 1j*randn(size(txSig)));
// plot(real(rxSig_demo), imag(rxSig_demo), 'b.');
// hold on;
// plot(real(txSig(1:100)), imag(txSig(1:100)), 'r+', 'LineWidth', 2); % Points idéaux
// title(['Constellation reçue à SNR = ' num2str(SNR_demo) ' dB']);
// xlabel('I'); ylabel('Q');
// legend('Reçu (Bruité)', 'Emis (Idéal)');
// grid on; axis square;


// %% Projet Télécom - Tâches I & II : Simulation Complète (Filtrée & Bruitée)
// clear; clc; close all;

// %% 1. Paramètres Globaux
// numBits = 50000;          % Nombre de bits (suffisant pour voir BER ~10^-4)
// SNR_dB_range = 0:2:16;    % Plage de SNR à tester
// sps = 4;                  % Échantillons par symbole (Oversampling)
// rolloff = 0.5;            % Facteur de retombée du filtre
// span = 6;                 % Longueur du filtre (en symboles)

// % Initialisation des résultats
// ber_qpsk = zeros(size(SNR_dB_range));
// ber_16qam = zeros(size(SNR_dB_range));

// % Création du filtre RRC (Root Raised Cosine)
// h = rcosdesign(rolloff, span, sps, 'sqrt');
// delay = span * sps / 2;   % Délai introduit par le filtre

// fprintf('--- Démarrage de la simulation comparative ---\n');

// %% 2. Boucle de Simulation sur le SNR
// for i = 1:length(SNR_dB_range)
//     SNR_dB = SNR_dB_range(i);
//     SNR_lin = 10^(SNR_dB/10);
    
//     %% --- CHAÎNE A : QPSK ---
//     M = 4; k = 2;
//     nBits_Q = ceil(numBits/k)*k;
//     dataIn_Q = randi([0 1], nBits_Q, 1);
    
//     % Modulation & Mapping
//     dataReshaped = reshape(dataIn_Q, k, length(dataIn_Q)/k)';
//     sym_I = 2*dataReshaped(:,1) - 1;
//     sym_Q = 2*dataReshaped(:,2) - 1;
//     txSym_Q = (sym_I + 1j*sym_Q) / sqrt(2);
    
//     % Sur-échantillonnage + Filtrage Tx
//     txUp_Q = upsample(txSym_Q, sps);
//     txSig_Q = conv(txUp_Q, h, 'same');
    
//     % Canal AWGN
//     % Calcul puissance signal filtré pour ajuster le bruit correctement
//     sigPower_Q = mean(abs(txSig_Q).^2);
//     noisePower_Q = sigPower_Q / SNR_lin;
//     noise_Q = sqrt(noisePower_Q/2) * (randn(size(txSig_Q)) + 1j*randn(size(txSig_Q)));
//     rxSig_Q = txSig_Q + noise_Q;
    
//     % Récepteur : Filtrage Rx + Sous-échantillonnage
//     rxFilt_Q = conv(rxSig_Q, h, 'same');
//     rxSym_Q = rxFilt_Q(1:sps:end); % Downsample
    
//     % Normalisation (Compensation gain filtre)
//     % Le filtre RRC en Tx + Rx change l'amplitude, on renormalise
//     rxSym_Q = rxSym_Q / sum(h.^2); 
    
//     % Démodulation
//     rx_bits_I = real(rxSym_Q) > 0;
//     rx_bits_Q = imag(rxSym_Q) > 0;
//     dataOut_Q = reshape([rx_bits_I rx_bits_Q]', nBits_Q, 1);
    
//     % BER QPSK
//     ber_qpsk(i) = sum(dataIn_Q ~= dataOut_Q) / nBits_Q;
    
//     %% --- CHAÎNE B : 16-QAM ---
//     M = 16; k = 4;
//     nBits_16 = ceil(numBits/k)*k;
//     dataIn_16 = randi([0 1], nBits_16, 1);
    
//     % Mapping 16-QAM (Manuel)
//     dataMat = reshape(dataIn_16, k, length(dataIn_16)/k)';
//     dec_I = dataMat(:,1)*2 + dataMat(:,2);
//     dec_Q = dataMat(:,3)*2 + dataMat(:,4);
//     lvl_I = 2*dec_I - 3; 
//     lvl_Q = 2*dec_Q - 3;
//     txSym_16 = (lvl_I + 1j*lvl_Q) / sqrt(10);
    
//     % Filtrage Tx
//     txUp_16 = upsample(txSym_16, sps);
//     txSig_16 = conv(txUp_16, h, 'same');
    
//     % Canal AWGN
//     sigPower_16 = mean(abs(txSig_16).^2);
//     noisePower_16 = sigPower_16 / SNR_lin;
//     noise_16 = sqrt(noisePower_16/2) * (randn(size(txSig_16)) + 1j*randn(size(txSig_16)));
//     rxSig_16 = txSig_16 + noise_16;
    
//     % Récepteur
//     rxFilt_16 = conv(rxSig_16, h, 'same');
//     rxSym_16 = rxFilt_16(1:sps:end);
//     rxSym_16 = rxSym_16 / sum(h.^2); % Normalisation filtre
//     rxSym_16 = rxSym_16 * sqrt(10);  % Dénormalisation énergie 16QAM
    
//     % Démodulation 16-QAM (Décision rigide)
//     rx_I_raw = real(rxSym_16); rx_Q_raw = imag(rxSym_16);
//     dec_I_RX = zeros(size(rx_I_raw)); dec_Q_RX = zeros(size(rx_Q_raw));
    
//     % Seuils : -2, 0, 2
//     dec_I_RX(rx_I_raw < -2) = 0;
//     dec_I_RX(rx_I_raw >= -2 & rx_I_raw < 0) = 1;
//     dec_I_RX(rx_I_raw >= 0 & rx_I_raw < 2) = 2;
//     dec_I_RX(rx_I_raw >= 2) = 3;
    
//     dec_Q_RX(rx_Q_raw < -2) = 0;
//     dec_Q_RX(rx_Q_raw >= -2 & rx_Q_raw < 0) = 1;
//     dec_Q_RX(rx_Q_raw >= 0 & rx_Q_raw < 2) = 2;
//     dec_Q_RX(rx_Q_raw >= 2) = 3;
    
//     b1_I = floor(dec_I_RX/2); b2_I = mod(dec_I_RX, 2);
//     b1_Q = floor(dec_Q_RX/2); b2_Q = mod(dec_Q_RX, 2);
//     dataOut_16 = reshape([b1_I b2_I b1_Q b2_Q]', nBits_16, 1);
    
//     % BER 16-QAM
//     ber_16qam(i) = sum(dataIn_16 ~= dataOut_16) / nBits_16;
// end

// %% 3. Visualisation Finale
// figure('Name', 'Performance BER avec Filtrage RRC');
// semilogy(SNR_dB_range, ber_qpsk, 'b-o', 'LineWidth', 2, 'DisplayName', 'QPSK (Filtré)');
// hold on;
// semilogy(SNR_dB_range, ber_16qam, 'r-s', 'LineWidth', 2, 'DisplayName', '16-QAM (Filtré)');
// grid on;
// xlabel('SNR (dB)');
// ylabel('BER (Bit Error Rate)');
// title('Comparaison de robustesse : QPSK vs 16-QAM (Canal Filtré)');
// legend;
// ylim([1e-5 1]);

// fprintf('Simulation terminée.\n');
