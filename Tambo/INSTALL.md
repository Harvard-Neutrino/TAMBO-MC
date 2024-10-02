# INSTALL Instructions for PROPOSAL, taurunner, and Configuring with PyCall

1. Need to make a Python venv
2. `pip install proposal`
3. `git clone taurunner`
4. Run `pip install .` inside taurunner home directory
5. Start Julia
6. Activate TAMBO env
7. Set `env["PYTHON"] = which python` inside Python venv
8. Set `env["PYTHONHOME"] = ""` (not always necessary)
9. Build PyCall
10. Success?
