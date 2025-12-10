%% code ameliore avec affichage dans le terminal 


%% Projet Télécom - Analyse Complète : BER & Constellations
clear; clc; close all;

%% 1. Paramètres de Simulation
numBits = 100000;         % Nombre de bits (élevé pour précision)
SNR_dB_range = 0:2:12;    % Plage de SNR
sps = 4;                  % Sur-échantillonnage (Samples per Symbol)
rolloff = 0.5;            % Filtre RRC
span = 6;
h = rcosdesign(rolloff, span, sps, 'sqrt'); % Filtre de mise en forme

% Variables pour stocker les résultats
ber_qpsk = zeros(size(SNR_dB_range));
ber_16qam = zeros(size(SNR_dB_range));

% Préparation de l'affichage Terminal
fprintf('---------------------------------------------------\n');
fprintf(' SNR (dB) |  BER QPSK  |  BER 16-QAM |  Erreurs  \n');
fprintf('---------------------------------------------------\n');

%% 2. Boucle de Simulation
for i = 1:length(SNR_dB_range)
    SNR_dB = SNR_dB_range(i);
    SNR_lin = 10^(SNR_dB/10);
    
    % --- QPSK ---
    M=4; k=2; nBits=ceil(numBits/k)*k;
    dataIn = randi([0 1], nBits, 1);
    % Modulation & Filtrage
    sym = (2*reshape(dataIn,k,[])' - (1+1i)) / sqrt(2); % QPSK rapide
    tx = conv(upsample(sym(:,1)+1j*sym(:,2), sps), h, 'same');
    % Canal
    noise = sqrt(mean(abs(tx).^2)/SNR_lin/2) * (randn(size(tx))+1j*randn(size(tx)));
    rx = conv(tx+noise, h, 'same');
    rx = rx(1:sps:end) / sum(h.^2);
    % BER
    dataOut = reshape([real(rx)>0 imag(rx)>0]', [], 1);
    err_qpsk = sum(dataIn~=dataOut);
    ber_qpsk(i) = err_qpsk / nBits;
    
    % Sauvegarde pour plot constellation (à 12 dB par exemple)
    if SNR_dB == 12
        rx_qpsk_plot = rx;
    end
    
    % --- 16-QAM ---
    M=16; k=4; nBits=ceil(numBits/k)*k;
    dataIn = randi([0 1], nBits, 1);
    % Mapping manuel rapide (Gray like)
    % Logique: conversion binaire -> niveaux PAM (-3,-1,1,3)
    bits = reshape(dataIn, k, [])';
    lvl_I = 2*(bits(:,1)*2+bits(:,2)) - 3;
    lvl_Q = 2*(bits(:,3)*2+bits(:,4)) - 3;
    sym = (lvl_I + 1j*lvl_Q) / sqrt(10);
    
    tx = conv(upsample(sym, sps), h, 'same');
    noise = sqrt(mean(abs(tx).^2)/SNR_lin/2) * (randn(size(tx))+1j*randn(size(tx)));
    rx = conv(tx+noise, h, 'same');
    rx = rx(1:sps:end) / sum(h.^2);
    rx_denorm = rx * sqrt(10); % Dénormalisation
    
    % Décision Rigide
    dec_I = zeros(size(rx)); dec_Q = zeros(size(rx));
    rI=real(rx_denorm); rQ=imag(rx_denorm);
    
    dec_I(rI<-2)=0; dec_I(rI>=-2 & rI<0)=1; dec_I(rI>=0 & rI<2)=2; dec_I(rI>=2)=3;
    dec_Q(rQ<-2)=0; dec_Q(rQ>=-2 & rQ<0)=1; dec_Q(rQ>=0 & rQ<2)=2; dec_Q(rQ>=2)=3;
    
    % Reconstitution bits
    bI1=floor(dec_I/2); bI2=mod(dec_I,2);
    bQ1=floor(dec_Q/2); bQ2=mod(dec_Q,2);
    dataOut = reshape([bI1 bI2 bQ1 bQ2]', [], 1);
    
    err_16qam = sum(dataIn~=dataOut);
    ber_16qam(i) = err_16qam / nBits;

    % Sauvegarde pour plot constellation
    if SNR_dB == 12
        rx_16qam_plot = rx;
    end
    
    % Affichage ligne par ligne
    fprintf('   %2d dB   |  %.2e  |  %.2e   | %d / %d \n', ...
        SNR_dB, ber_qpsk(i), ber_16qam(i), err_16qam, err_qpsk);
end
fprintf('---------------------------------------------------\n');

%% 3. Visualisation : Courbes BER
figure('Name', 'Résultats Tâche II', 'Position', [100, 100, 1000, 500]);

subplot(1, 2, 1);
semilogy(SNR_dB_range, ber_qpsk, 'b-o', 'LineWidth', 2, 'DisplayName', 'QPSK');
hold on;
semilogy(SNR_dB_range, ber_16qam, 'r-s', 'LineWidth', 2, 'DisplayName', '16-QAM');
grid on;
xlabel('SNR (dB)'); ylabel('BER');
title('Performance BER (Bit Error Rate)');
legend('Location', 'SouthWest');
ylim([1e-5 1]);

%% 4. Visualisation : Constellations Comparées (à SNR = 12 dB)
subplot(1, 2, 2);
plot(real(rx_16qam_plot), imag(rx_16qam_plot), 'r.', 'MarkerSize', 4);
hold on;
plot(real(rx_qpsk_plot), imag(rx_qpsk_plot), 'b.', 'MarkerSize', 4);
grid on; axis square;
title('Constellations Reçues @ SNR = 12dB');
legend('16-QAM (Bruité)', 'QPSK (Clair)');
xlabel('I'); ylabel('Q');

fprintf('Simulation terminée. Les figures sont prêtes.\n');