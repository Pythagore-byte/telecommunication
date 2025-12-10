%% Projet Télécom - Tâche III : WebLab avec Signal QPSK/16-QAM (Version Finale)
% Ce script génère un signal, l'envoie au banc d'essai à Barcelone,
% et analyse les résultats (EVM, ACPR, AM/AM) sans faire planter Matlab.

clear all; close all; clc;

%% 1. Configuration WebLab (Ne pas toucher)
Window='Blackman-Harris';
SegLength=2^10;
overlap=25;
Hs = spectrum.welch(Window,SegLength,overlap);
EnAMbL='on'; EnSpec='on';

%% 2. Génération de NOTRE Stimulus (Remplacement de 'GetStimulus')
fprintf('--- Génération du signal QPSK/16-QAM ---\n');

% --- Paramètres du signal ---
M = 16;                 % 4 = QPSK, 16 = 16-QAM
Fs = 200e6;             % Fréquence d'échantillonnage WebLab (200 MHz)
Bw = 10e6;              % Bande passante (10 MHz)
sps = 4;                % Échantillons par symbole

% Calculs dérivés
numSymb = 20000;        % Nombre de symboles
numBits = numSymb * log2(M);

% Génération des bits aléatoires
dataIn = randi([0 1], numBits, 1);

% --- Modulation (QPSK ou 16-QAM) ---
k = log2(M);
bits_reshaped = reshape(dataIn, k, [])';

if M == 4 % QPSK
    sym = (2*bits_reshaped(:,1)-1 + 1j*(2*bits_reshaped(:,2)-1))/sqrt(2);
elseif M == 16 % 16-QAM
    lvl_I = 2*(bits_reshaped(:,1)*2+bits_reshaped(:,2)) - 3;
    lvl_Q = 2*(bits_reshaped(:,3)*2+bits_reshaped(:,4)) - 3;
    sym = (lvl_I + 1j*lvl_Q)/sqrt(10);
end

% --- Filtrage RRC (Mise en forme) ---
rolloff = 0.5; span = 6;
h = rcosdesign(rolloff, span, sps, 'sqrt');
txSignal = conv(upsample(sym, sps), h, 'same');

% --- ADAPTATION WEBLAB ---
% 1. Normalisation (Amplitude max à 1)
PAin = txSignal / max(abs(txSignal)); 
% 2. Format Vecteur Colonne
PAin = PAin(:); 

% 3. Configuration ACPR pour les fonctions du prof
ACPR.BW = Bw;           
ACPR.Fs = Fs;           
ACPR.Offset = Bw;       % Décalage pour mesurer le canal adjacent (Important !)
ACPR.Name = 'Default';  

fprintf('Signal généré : %d échantillons.\n', length(PAin));

%% 3. Mesure sur l'Amplificateur Réel (Connexion WebLab)
fprintf('--- Envoi vers WebLab en cours... ---\n');

% Puissance Cible (RMSin)
% -22 dBm : Régime linéaire (Signal propre)
% -15 dBm : Régime saturé (Signal distordu -> Utile pour tester la DPD)
RMSin_Target = -14; 

try
    [PAout, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_2(PAin, RMSin_Target); 
catch
    error('Erreur de connexion. Vérifie ta connexion internet ou le VPN.');
end

% Protection si le WebLab rejette le signal
if isempty(PAout)
    error('Le WebLab a rejeté le signal (Power check failed). Baisse RMSin_Target.');
end

fprintf('Réception terminée.\n');
fprintf('Puissance Sortie Mesurée : %.2f dBm\n', RMSout);

%% 4. Analyse et Visualisation
% A. Alignement Temporel (Correction du délai de transmission)
PAout = timealign(PAin, PAout);

% B. Tracé AM/AM et AM/PM (Optimisé)
Data.In = PAin;
Data.Out = PAout;

Data.ALimLinIn = 0.2; 

% --- OPTIMISATION ANTI-CRASH ---
% On ne garde qu'un point sur 50 pour l'affichage graphique
% (Sinon Matlab sature et affiche des warnings rouges)
step_plot = 50; 
Data.In = Data.In(1:step_plot:end);
Data.Out = Data.Out(1:step_plot:end);

figure('Name', 'AM/AM & AM/PM');
PA_Metrics = amampm(Data, EnAMbL); % Un seul appel ici !

%% C. Analyse Spectrale (ACPR) - VERSION CORRIGÉE
figure('Name', 'Spectres Entrée/Sortie');

% 1. Calculs des spectres
[ACPRin, PSDin] = acpr(PAin, Fs, ACPR);    
[ACPRout, PSDout] = acpr(PAout, Fs, ACPR); 

% 2. Extraction des données pour le tracé manuel
% (On s'assure de récupérer les vecteurs de fréquence et de puissance)
f_axis = PSDin.Frequencies / 1e6; % Conversion en MHz pour l'affichage
P_in_dB = 10*log10(PSDin.Data);   % Conversion en dB
P_out_dB = 10*log10(PSDout.Data);

% 3. Tracé Manuel (Plot)
plot(f_axis, P_in_dB, 'b', 'LineWidth', 1.5, 'DisplayName', 'Entrée (PAin)');
hold on;
plot(f_axis, P_out_dB, 'r', 'LineWidth', 1.5, 'DisplayName', 'Sortie (PAout)');

% 4. Mise en forme du graphique
grid on;
xlabel('Fréquence (MHz)');
ylabel('Densité Spectrale de Puissance (dB)');
legend show;
title(['Spectres @ Pin = ' num2str(RMSin_Target) ' dBm']);
ylim([-100 0]); % Zoom sur la partie intéressante (ajuste si nécessaire)
%% 5. Calcul des Métriques (EVM & ACPR)
% Recalcul du gain complexe sur le signal complet (pas celui réduit pour le plot)
gain_complexe = (PAin' * PAout) / (PAin' * PAin);
PAout_norm = PAout / gain_complexe;

% Calcul EVM
error_vec = abs(PAin - PAout_norm).^2;
evm_rms = sqrt(mean(error_vec) / mean(abs(PAin).^2)) * 100;

fprintf('\n--- Résultats Finaux ---\n');
fprintf('EVM (Error Vector Magnitude) : %.2f %%\n', evm_rms);
fprintf('ACPR (Canal Supérieur)       : %.2f dB\n', ACPRout.U1);
fprintf('ACPR (Canal Inférieur)       : %.2f dB\n', ACPRout.L1);