%% Projet Télécom - Tâche III : WebLab & DPD (Version Professeur)
clear all; close all; clc;

Window='Blackman-Harris';
SegLength=2^10;
overlap=25;
Hs = spectrum.welch(Window,SegLength,overlap);
EnAMbL='on'; EnSpec='on';
    

%% 1. Configuration & Génération Signal
% On garde notre générateur qui marche bien
M = 16; Fs = 200e6; Bw = 10e6; sps = 4;
RMSin_Target = -15; % Puissance de saturation (ajuste si besoin)

% Génération QPSK/16QAM
fprintf('--- 1. Génération du signal ---\n');
numSymb = 20000; numBits = numSymb * log2(M);
dataIn = randi([0 1], numBits, 1);

k = log2(M); bits_reshaped = reshape(dataIn, k, [])';
if M==16 
    lvl_I = 2*(bits_reshaped(:,1)*2+bits_reshaped(:,2))-3; 
    lvl_Q = 2*(bits_reshaped(:,3)*2+bits_reshaped(:,4))-3;
    sym = (lvl_I+1j*lvl_Q)/sqrt(10);
else
    sym = (2*bits_reshaped(:,1)-1 + 1j*(2*bits_reshaped(:,2)-1))/sqrt(2);
end

% Filtrage et Normalisation
rolloff = 0.5; span = 6;
h = rcosdesign(rolloff, span, sps, 'sqrt');
txSignal = conv(upsample(sym, sps), h, 'same');
PAin = txSignal / max(abs(txSignal)); % Normalisation crête à 1
PAin = PAin(:);

% Configuration ACPR (Requise par la fonction du prof)
ACPR.BW = Bw; ACPR.Fs = Fs; ACPR.Offset = Bw; ACPR.Name = 'Default';

%% 2. Mesure SANS DPD (Référence)
fprintf('--- 2. Envoi SANS DPD (Référence) ---\n');
try
    [PAout, RMSout] = RFWebLab_PA_meas_v1_2(PAin, RMSin_Target);
    PAout = timealign(PAin, PAout);
catch
    error('Erreur WebLab. Vérifie ta connexion/VPN.');
end

% Calcul ACPR initial (pour donner à la fonction DPD)
[ACPRout, PSDout] = acpr(PAout, Fs, ACPR);

%% 3. Configuration pour la fonction DPD du Prof
fprintf('--- 3. Lancement de l''algorithme DPD (Fichier Prof) ---\n');

% A. Structure 'Algorithm'
Algorithm.Method = 'Alternate'; % Ou 'Backward', selon ce qui marche le mieux
Algorithm.NbIterPerStage = 3;   % Nombre d'itérations (3 suffisent souvent)
Algorithm.NbSampPerIter = length(PAin); % On utilise tout le signal pour apprendre
Algorithm.DampingNewtonFactor = 0.7; % Facteur d'apprentissage (mu)

% B. Structure 'PD' (Données d'entrée)
PD.In = PAin;
PD.BW = Bw;

% C. Structure 'Stimulus'
Stimulus.wf = PAin;
Stimulus.RMSin = RMSin_Target;
Stimulus.Fs = Fs;
Stimulus.ACPR = ACPR;

%% 4. Appel de la fonction DPD.m
% C'est ici qu'on utilise le fichier que tu m'as donné
% La fonction va ouvrir ses propres figures et tourner en boucle
SystIter = DPD(Stimulus, ACPRout, Algorithm, PD);

fprintf('DPD terminée.\n');

%% 5. Analyse des Résultats Finaux
% La fonction DPD renvoie une structure 'SystIter'. 
% La dernière case contient le meilleur résultat.

LastIter = SystIter(end);       % Dernière itération
PAout_DPD = LastIter.y;         % Signal de sortie linéarisé
PAin_DPD = LastIter.PD.x;       % Signal pré-distordu envoyé

% Alignement final pour calculs
PAout_DPD = timealign(PAin, PAout_DPD);

% --- Calcul EVM ---
gain_ref = (PAin'*PAout)/(PAin'*PAin);
evm_no_dpd = sqrt(mean(abs(PAin - PAout/gain_ref).^2)/mean(abs(PAin).^2))*100;

gain_dpd = (PAin'*PAout_DPD)/(PAin'*PAin);
evm_dpd = sqrt(mean(abs(PAin - PAout_DPD/gain_dpd).^2)/mean(abs(PAin).^2))*100;

fprintf('\n--- RÉSULTATS FINAUX ---\n');
fprintf('EVM Sans DPD : %.2f %%\n', evm_no_dpd);
fprintf('EVM Avec DPD : %.2f %%\n', evm_dpd);
fprintf('Gain ACPR    : %.2f dB\n', LastIter.ACPRimpr.U1); % Amélioration canal adjacent

%% 6. Comparaison Graphique (AM/AM)
figure('Name', 'Comparaison Finale AM/AM');
step = 50; % Optimisation affichage
plot(abs(PAin(1:step:end)), abs(PAout(1:step:end)/gain_ref), 'r.', 'DisplayName', 'Sans DPD');
hold on;
plot(abs(PAin(1:step:end)), abs(PAout_DPD(1:step:end)/gain_dpd), 'g.', 'DisplayName', 'Avec DPD');
plot([0 1], [0 1], 'k--', 'LineWidth', 2);
legend; grid on; title('AM/AM : Avant et Après DPD'); xlabel('Entrée'); ylabel('Sortie');

%% 7. Comparaison Spectrale
figure('Name', 'Comparaison Finale Spectres');
f_axis = PSDout.Frequencies/1e6;
plot(f_axis, 10*log10(PSDout.Data), 'r', 'LineWidth', 1.5, 'DisplayName', 'Sans DPD');
hold on;
% On récupère le spectre de la dernière itération
PSD_DPD_Data = 10.^(LastIter.PY/10); 
plot(f_axis, 10*log10(PSD_DPD_Data), 'g', 'LineWidth', 1.5, 'DisplayName', 'Avec DPD');
grid on; legend; title('Spectres : Réduction de la repousse'); xlabel('MHz'); ylabel('dB'); ylim([-80 0]);