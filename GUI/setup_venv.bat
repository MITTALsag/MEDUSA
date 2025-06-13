@echo off

REM Nom de l'environnement virtuel
set ENV_NAME=pyqt6_venv

REM Chemin où l'environnement virtuel sera créé
set ENV_DIR=.\%ENV_NAME%

REM Créer l'environnement virtuel
python3 -m venv %ENV_DIR%

REM Activer l'environnement virtuel
call %ENV_DIR%\Scripts\activate.bat

REM Installer PyQt6
pip3 install PyQt6

echo Environnement virtuel créé et PyQt6 installé dans %ENV_DIR%

REM Désactiver l'environnement virtuel
deactivate

call activer_venv.bat