%% Projet Télécom - Étape 1 : Chaîne de transmission de base (Sans Bruit)
clear; clc; close all;

%% 1. Paramètres de simulation
M = 4;                  % Ordre de modulation (4 = QPSK, 16 = 16QAM) [cite: 76]
k = log2(M);            % Nombre de bits par symbole
numBits = 10000;        % Nombre total de bits à transmettre (doit être multiple de k)
numBits = ceil(numBits/k)*k; % Ajustement pour être multiple de k

%% 2. Émetteur (Transmitter TX) [cite: 74]
% Génération des données binaires aléatoires
dataIn = randi([0 1], numBits, 1); 

% Modulation (Mappage binaire vers symboles complexes)
% On utilise 'UnitAveragePower' pour que la puissance moyenne soit de 1 Watt
dataMod = qammod(dataIn, M, 'InputType', 'bit', 'UnitAveragePower', true);

%% 3. Canal (Idéal pour l'instant) [cite: 71]
% Pas de bruit, pas de délai, pas d'atténuation
receivedSignal = dataMod;

%% 4. Récepteur (Receiver RX) [cite: 80]
% Démodulation
dataOut = qamdemod(receivedSignal, M, 'OutputType', 'bit', 'UnitAveragePower', true);

%% 5. Analyse et Visualisation [cite: 83]

% a) Calcul du BER (Bit Error Rate) 
[number_of_errors, bit_error_rate] = biterr(dataIn, dataOut);

fprintf('--- Résultats de la simulation ---\n');
fprintf('Modulation : %d-QAM\n', M);
fprintf('Nombre d''erreurs : %d\n', number_of_errors);
fprintf('Taux d''erreur binaire (BER) : %.4e\n', bit_error_rate);

% b) Visualisation de la constellation [cite: 86]
figure;
scatterplot(receivedSignal);
title(['Constellation reçue (Sans bruit) - ' num2str(M) '-QAM']);
grid on;
