%% Projet Télécom - Tâche III (Partie 1) : Caractérisation de l'Ampli
% Ce script mesure tous les paramètres physiques demandés avant la DPD.
clear all; close all; clc;

%% 1. Configuration
Window='Blackman-Harris'; SegLength=2^10; overlap=25;
Hs = spectrum.welch(Window,SegLength,overlap);
EnAMbL='on'; EnSpec='on';

% Paramètres
M = 16; Fs = 200e6; Bw = 10e6; sps = 4;

% --- CHOIX DE LA PUISSANCE ---
% Change cette valeur pour faire tes deux tests :
% Test 1 (Linéaire) : -22
% Test 2 (Saturé)   : -15
RMSin_Target = -14; 

%% 2. Génération du Signal (16-QAM)
fprintf('--- Génération 16-QAM ---\n');
numSymb = 20000; numBits = numSymb * log2(M);
dataIn = randi([0 1], numBits, 1);

k = log2(M); bits_reshaped = reshape(dataIn, k, [])';
lvl_I = 2*(bits_reshaped(:,1)*2+bits_reshaped(:,2))-3; 
lvl_Q = 2*(bits_reshaped(:,3)*2+bits_reshaped(:,4))-3;
sym = (lvl_I+1j*lvl_Q)/sqrt(10);

h = rcosdesign(0.5, 6, sps, 'sqrt');
txSignal = conv(upsample(sym, sps), h, 'same');
PAin = txSignal / max(abs(txSignal)); % Normalisation
PAin = PAin(:);

ACPR.BW = Bw; ACPR.Fs = Fs; ACPR.Offset = Bw; ACPR.Name = 'Default';

%% 3. Mesure WebLab (Envoi du signal)
fprintf('--- Envoi au WebLab (@ %.1f dBm) ---\n', RMSin_Target);
try
    % La fonction renvoie aussi le Courant (Idc) et la Tension (Vdc) !
    [PAout, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_2(PAin, RMSin_Target);
    PAout = timealign(PAin, PAout);
catch
    error('Erreur connexion WebLab.');
end

%% 4. Calcul des Métriques Électriques & Efficacité
% Formule de l'efficacité ajoutée (PAE - Power Added Efficiency)
% P_out (Watts) = 10^(dBm/10) / 1000
P_out_Watt = 10^(RMSout/10) / 1000;
P_in_Watt  = 10^(RMSin_Target/10) / 1000;
P_DC_Watt  = Vdc * Idc; % Puissance consommée par l'alim

% Rendement (Drain Efficiency) = P_RF_out / P_DC
Efficiency = (P_out_Watt / P_DC_Watt) * 100;

fprintf('\n--- Métriques Électriques ---\n');
fprintf('Tension Alim (Vdc) : %.2f V\n', Vdc);
fprintf('Courant Alim (Idc) : %.2f A\n', Idc);
fprintf('Puissance DC       : %.2f W\n', P_DC_Watt);
fprintf('Puissance RF Out   : %.2f W (%.2f dBm)\n', P_out_Watt, RMSout);
fprintf('EFFICACITÉ (Drain) : %.2f %%\n', Efficiency);

%% 5. Calcul EVM et ACPR
gain_comp = (PAin'*PAout)/(PAin'*PAin);
PAout_norm = PAout / gain_comp;
evm = sqrt(mean(abs(PAin - PAout_norm).^2)/mean(abs(PAin).^2))*100;
[ACPRout, ~] = acpr(PAout, Fs, ACPR);

fprintf('\n--- Qualité Signal ---\n');
fprintf('EVM  : %.2f %%\n', evm);
fprintf('ACPR : %.2f dB (Haut) / %.2f dB (Bas)\n', ACPRout.U1, ACPRout.L1);

%% 6. Calcul BER (Avec Récepteur complet)
rxFilt = conv(PAout, h, 'same');
rxSym = rxFilt(1:sps:end);
% Normalisation AGC précise
rxSym = rxSym / (sum(h.^2) * gain_comp);

% Démod 16-QAM
rxSym = rxSym * sqrt(10);
rI=real(rxSym); rQ=imag(rxSym);
dI=zeros(size(rI)); dQ=zeros(size(rQ));
dI(rI<-2)=0; dI(rI>=-2&rI<0)=1; dI(rI>=0&rI<2)=2; dI(rI>=2)=3;
dQ(rQ<-2)=0; dQ(rQ>=-2&rQ<0)=1; dQ(rQ>=0&rQ<2)=2; dQ(rQ>=2)=3;
bI1=floor(dI/2); bI2=mod(dI,2); bQ1=floor(dQ/2); bQ2=mod(dQ,2);
rx_bits = reshape([bI1 bI2 bQ1 bQ2]', [], 1);

err = sum(dataIn(:) ~= rx_bits(:));
fprintf('BER  : %.2e (%d erreurs)\n', err/length(dataIn), err);

%% 7. Visualisations Complètes
% A. AM/AM
Data.In = PAin(1:50:end); Data.Out = PAout(1:50:end); Data.ALimLinIn=0.2;
figure('Name', 'AM/AM'); amampm(Data, EnAMbL);

% B. Constellation Dégradée (Ce que demande le PDF)
figure('Name', 'Constellation Reçue');
plot(real(rxSym(1:2000)), imag(rxSym(1:2000)), 'r.');
hold on;
% Points idéaux (centres)
ref_pts = [-3 -1 1 3];
[GX,GY] = meshgrid(ref_pts, ref_pts);
plot(GX(:), GY(:), 'k+', 'MarkerSize', 10, 'LineWidth', 2);
title(['Constellation Reçue @ ' num2str(RMSin_Target) ' dBm']);
grid on; axis square;
xlabel('I'); ylabel('Q');