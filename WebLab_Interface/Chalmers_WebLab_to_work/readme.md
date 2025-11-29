# Projet Télécom : Simulateur de Transmission & Hardware-in-the-Loop (WebLab)

**Université :** Sorbonne Université  
**Cours :** Projet EI Télécom S9  
**Technologies :** Matlab, Signal Processing Toolbox, RF WebLab (Hardware-in-the-loop)

---

## 📋 Description du Projet
Ce projet vise à concevoir une chaîne de transmission numérique complète (émetteur-récepteur) et à étudier les effets des non-linéarités d'un Amplificateur de Puissance (PA) réel.

Le projet se divise en deux parties principales :
1.  **Simulation Matlab :** Implémentation d'une chaîne TX/RX avec modulations QPSK/16-QAM, filtrage de mise en forme et canal bruité (AWGN).
2.  **Mesure Réelle (WebLab) :** Pilotage à distance d'un banc de test RF (situé à Barcelone) pour caractériser un amplificateur de puissance réel et observer ses défauts (compression, repousse spectrale).

---

## 🚀 Fonctionnalités Implémentées

### 1. Simulation de la Chaîne de Transmission (Matlab)
* **Modulations :** QPSK et 16-QAM.
* **Mise en forme :** Filtre en Racine de Cosinus Surélevé (Root Raised Cosine - RRC) pour limiter la bande passante et minimiser l'IES.
* **Canal :** Ajout de bruit blanc gaussien (AWGN).
* **Analyse de Performance :**
    * Calcul du **BER** (Bit Error Rate) en fonction du SNR.
    * Comparaison théorique vs simulée.
    * Visualisation des constellations (diagrammes de l'œil).
* **Démonstrateurs :**
    * Transmission de texte ASCII.
    * Transmission d'image (Bitmap) à travers le canal bruité.

### 2. Caractérisation RF (WebLab)
* Connexion au serveur distant WebLab (UPC).
* Génération et envoi de signaux QPSK/16-QAM formatés pour le banc de test.
* **Mesures effectuées :**
    * **EVM** (Error Vector Magnitude) : Mesure de la qualité du signal.
    * **ACPR** (Adjacent Channel Power Ratio) : Mesure de la pollution spectrale.
    * **Courbes AM/AM et AM/PM** : Caractérisation de la non-linéarité et de la saturation de l'amplificateur.

---

## 📂 Structure des Fichiers

* `projet_telecom_final_step2.m` : Script principal de simulation Matlab. Génère les courbes de BER comparatives (QPSK vs 16-QAM) et valide la chaîne théorique.
* `projet_telecom_demo_image.m` : Démonstrateur ludique transmettant une image pixel par pixel à travers le canal bruité QPSK.
* `main_weblab_final.m` : Script d'interface avec le WebLab. Il génère le signal, l'envoie à l'ampli réel, récupère la sortie et trace les caractéristiques AM/AM et spectrales.
* `weblab_files/` : Contient les fonctions fournies pour la connexion (`RFWebLab_PA_meas_v1_2.m`, `timealign.m`, etc.).

---

## 📊 Résultats Obtenus (Aperçu)

### 1. Robustesse au Bruit
Nous avons validé que la modulation **QPSK** est plus robuste que la **16-QAM** face au bruit.
* *Seuil de performance :* La QPSK atteint un BER de $10^{-4}$ autour de 6 dB de SNR, contre 13 dB pour la 16-QAM.

### 2. Impact du Bruit sur l'Image
Visualisation de l'effet du canal AWGN sur une image transmise en QPSK :
* **SNR Fort (12 dB) :** Image parfaite, constellation nette.
* **SNR Faible (3 dB) :** Image bruitée ("poivre et sel"), constellation dispersée.

### 3. Non-linéarité de l'Amplificateur (WebLab)
En poussant la puissance d'entrée à **-14 dBm**, nous avons mis en évidence la saturation de l'amplificateur réel :
* **Compression de Gain :** Visible sur la courbe AM/AM (aplatissement aux fortes amplitudes).
* **Repousse Spectrale :** Apparition d'épaules sur le spectre de sortie (ACPR dégradé).
* **EVM :** Augmentation significative de l'erreur (de ~4% en régime linéaire à >8% en saturation).

---

## 🛠️ Prochaines Étapes
* Implémentation de l'algorithme de **Linéarisation par Prédistorsion Numérique (DPD)**.
* Utilisation de modèles polynômiaux à mémoire (Memory Polynomials).
* Validation de la correction sur le banc WebLab.

---
**Auteurs :** [TOURE SEKOUBA] & [ARTHUR ] & [KARLITOU]