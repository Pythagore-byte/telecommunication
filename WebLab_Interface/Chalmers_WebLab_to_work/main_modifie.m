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
    RMSin_Target = -15; 
    
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
    %% 6. Calcul du BER avec Distorsion PA (Tâche III)
    fprintf('\n--- Performance BER (Bit Error Rate) ---\n');
    
    % --- CORRECTION : AJOUT DU RÉCEPTEUR NUMÉRIQUE ---
    % PAout est le signal temporel (sur-échantillonné par 4).
    % Il faut le filtrer et le sous-échantillonner pour retrouver les symboles.
    
    % 1. Filtrage Adapté (Matched Filter)
    rxFiltered_PA = conv(PAout, h, 'same');
    
    % 2. Sous-échantillonnage (Downsampling) : On garde 1 point sur 4
    rxSym_PA = rxFiltered_PA(1:sps:end);
    
    % 3. Normalisation
    % On compense le gain du filtre RRC
    rxSym_PA = rxSym_PA / sum(h.^2);
    
    % On compense le gain de l'ampli (AGC) pour avoir une puissance unitaire
    % On utilise un calcul simple de gain moyen entre émis et reçu
    gain_est = (sym' * rxSym_PA) / (sym' * sym);
    rxSym_PA = rxSym_PA / gain_est;
    
    % --- DÉMODULATION ---
    rx_bits = zeros(length(dataIn), 1);
    
    if M == 4 % QPSK
        rx_bits_I = real(rxSym_PA) > 0;
        rx_bits_Q = imag(rxSym_PA) > 0;
        rx_bits = reshape([rx_bits_I rx_bits_Q]', [], 1);
        
    elseif M == 16 % 16-QAM
        % Dénormalisation pour retrouver les niveaux entiers (-3, -1, 1, 3)
        rx_denorm = rxSym_PA * sqrt(10); 
        
        rI = real(rx_denorm); rQ = imag(rx_denorm);
        dec_I = zeros(size(rI)); dec_Q = zeros(size(rQ));
        
        % Seuils de décision 16-QAM (-2, 0, 2)
        dec_I(rI<-2)=0; dec_I(rI>=-2 & rI<0)=1; dec_I(rI>=0 & rI<2)=2; dec_I(rI>=2)=3;
        dec_Q(rQ<-2)=0; dec_Q(rQ>=-2 & rQ<0)=1; dec_Q(rQ>=0 & rQ<2)=2; dec_Q(rQ>=2)=3;
        
        bI1=floor(dec_I/2); bI2=mod(dec_I,2);
        bQ1=floor(dec_Q/2); bQ2=mod(dec_Q,2);
        rx_bits = reshape([bI1 bI2 bQ1 bQ2]', [], 1);
    end
    
    % 4. Comparaison (Maintenant ils ont la même taille !)
    % On s'assure que les vecteurs sont bien colonnes
    rx_bits = rx_bits(:); 
    dataIn = dataIn(:);
    
    % Sécurité taille (au cas où un léger décalage subsiste)
    L = min(length(dataIn), length(rx_bits));
    numErrors = sum(dataIn(1:L) ~= rx_bits(1:L));
    BER_PA = numErrors / L;
    
    fprintf('Bits envoyés : %d\n', L);
    fprintf('Erreurs      : %d\n', numErrors);
    fprintf('BER Mesuré   : %.2e\n', BER_PA);
    
    if BER_PA == 0
        fprintf('=> Transmission parfaite (Zone Linéaire)\n');
    else
        fprintf('=> Dégradation due à la non-linéarité (Zone de Saturation)\n');
    end

    