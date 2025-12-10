%% Projet Télécom - Démo : Transmission d'Image sur canal Bruité
clear; clc; close all;

%% 1. Configuration et Chargement de l'Image
% --- Paramètres ---
SNR_dB = -10;  % Essaie 5 dB (très bruité), 8 dB (moyen), 12 dB (parfait)
% ------------------

% Chargement d'une image intégrée à Matlab 
try
    img_original = imread('messi.jpg'); % Image noir et blanc standard
    %  pour ma propre image, je peux utiliser la ligne dessous :
    % img_original = imread('mon_image.jpg'); 
catch
    % Si l'image n'existe pas, on crée un damier simple
    img_original = checkerboard(20) > 0.5; 
    img_original = uint8(img_original * 255);
end

% On redimensionne pour que la simulation soit rapide (ex: 100x100 pixels)
% Une image trop grande prendrait trop de temps à simuler bit par bit
scale_factor = 0.5; % Réduire l'image (0.5 = 50% de la taille)
img_resized = imresize(img_original, scale_factor); 
%cette fonctionne permet de redimensionner l'image
[rows, cols, channels] = size(img_resized);  % on cree 3 variables pour stocker les parametres de l'image a savoir nombre de colonne , de ligne , et nombre de  couleur disponible dans l'image

fprintf('Image chargée : %dx%d pixels\n', rows, cols);

%% 2. Conversion Image -> Bits (Source Binaire)
% Aplatir l'image en un vecteur de pixels
pixels_vec = img_resized(:);  % mettre tous les pixels de l image sur une seule colonne . ( empiler les uns apres les autres)

% Convertir en binaire (8 bits par pixel)
% dec2bin convertit en char, on repasse en numérique
bin_matrix = dec2bin(pixels_vec, 8); % convertir chaque pixel(0-255) en bit (les 0 et 1)
bits_temp = bin_matrix.' - '0'; % Transpose pour lire bit par bit (en colonne) , et on soustrait -'0' pour convertir le texte brut en nombre 
dataIn = bits_temp(:);

numBits = length(dataIn);
fprintf('Transmission de %d bits en cours...\n', numBits);

%% 3. Chaîne de Transmission (QPSK + Filtrage RRC)
% Paramètres
M = 4; k = 2; % QPSK
sps = 4; rolloff = 0.5; span = 6;
h = rcosdesign(rolloff, span, sps, 'sqrt'); 

% --- MODULATEUR ---
% Padding si nécessaire (multiple de k)
remainder = mod(numBits, k);
if remainder > 0
    padding = k - remainder;
    dataIn = [dataIn; zeros(padding, 1)];
end

% Modulation QPSK
bits_reshaped = reshape(dataIn, k, [])';
sym_I = 2*bits_reshaped(:,1) - 1;
sym_Q = 2*bits_reshaped(:,2) - 1;
sym = (sym_I + 1j*sym_Q) / sqrt(2);

% Filtrage émission
txSignal = conv(upsample(sym, sps), h, 'same');

% --- CANAL AWGN ---
SNR_lin = 10^(SNR_dB/10);
sigPower = mean(abs(txSignal).^2);
noisePower = sigPower / SNR_lin;
noise = sqrt(noisePower/2) * (randn(size(txSignal)) + 1j*randn(size(txSignal)));
rxSignal = txSignal + noise;

% --- RÉCEPTEUR ---
% Filtrage réception
rxFiltered = conv(rxSignal, h, 'same');
rxSym = rxFiltered(1:sps:end) / sum(h.^2);

% Démodulation (Décision rigide)
rx_bits_I = real(rxSym) > 0;
rx_bits_Q = imag(rxSym) > 0;
dataOut_temp = reshape([rx_bits_I rx_bits_Q]', [], 1);

% Retirer le padding
dataOut = dataOut_temp(1:numBits);

%% 4. Reconstruction de l'Image et Calcul BER
% Calcul du BER
nb_erreurs = sum(dataIn ~= dataOut);
BER = nb_erreurs / numBits;

% Conversion Bits -> Pixels
% On remet les bits sous forme de matrice (8 lignes)
bits_rx_mat = reshape(dataOut, 8, []).';
% Conversion binaire -> entier (0-255)
pixels_rx_vec = bin2dec(char(bits_rx_mat + '0'));

% Conversion vecteur -> Image (uint8)
img_reconstructed = reshape(uint8(pixels_rx_vec), rows, cols, channels);

fprintf('Réception terminée.\n');
fprintf('SNR : %d dB | BER : %.4f (%d erreurs)\n', SNR_dB, BER, nb_erreurs);

%% 5. Visualisation
figure('Name', 'Transmission Image QPSK', 'Position', [100, 100, 1000, 400]);

% Image Originale
subplot(1, 3, 1);
imshow(img_resized);
title('1. Image Envoyée (Originale)');

% Image Reçue
subplot(1, 3, 2);
imshow(img_reconstructed);
title({['2. Image Reçue (SNR = ' num2str(SNR_dB) ' dB)'], ...
       ['BER = ' num2str(BER, '%.4f')]});

% Constellation
subplot(1, 3, 3);
plot(real(rxSym), imag(rxSym), 'b.', 'MarkerSize', 1);
hold on;
plot(real(sym(1:100)), imag(sym(1:100)), 'r+', 'LineWidth', 2); % Quelques points idéaux
% Au lieu de prendre les 100 premiers, on prend les valeurs uniques idéales
symboles_uniques = unique(sym); 
plot(real(symboles_uniques), imag(symboles_uniques), 'r+', 'LineWidth', 2, 'MarkerSize', 15);

axis square; grid on;
title('3. Constellation Reçue');
xlabel('I'); ylabel('Q');