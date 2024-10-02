# INSTALL Instructions for PROPOSAL, taurunner, and Configuring with PyCall

1. Need to make a Python venv (`python3 -m venv your-env-name`, then `source your-env-name/activate/bin`)
2. `pip install proposal`
3. `git clone` taurunner 
4. Run `pip install .` inside `path/to/taurunner`
5. Start Julia
6. Activate TAMBO `path/to/tambo/tambo`
7. Set `env["PYTHON"]` = the path from `which python` inside Python venv
8. Set `env["PYTHONHOME"] = ""` and `env["PYTHONPATH"]=""` (not always necessary)
9. `Using Pkg; using Tambo; Pkg.build("PyCall")` Build PyCall
10. Success? type `using PyCall; tr = pyimport("taurunner")` to see if successful 
