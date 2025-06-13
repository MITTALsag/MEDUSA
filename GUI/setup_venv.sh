#!/bin/bash

# Nom de l'environnement virtuel
ENV_NAME="pyqt6_venv"

# Chemin où l'environnement virtuel sera créé
ENV_DIR="./$ENV_NAME"

# Créer l'environnement virtuel
python3 -m venv $ENV_DIR

# Activer l'environnement virtuel
source $ENV_DIR/bin/activate

# Installer PyQt6
pip3 install pyqt6

echo "Environnement virtuel créé et PyQt6 installé dans $ENV_DIR"

deactivate

./activer_venv.sh