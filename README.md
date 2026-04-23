This code is supposed to run a comparison between a symplectic integrator and RK4 that are both applied to the Jovian moons (currently the Galilean moons). It uses JPL's ephemerides as a source of truth.

## Requirements
- Python 3.14.x
- numpy
- requests 2.32.5
- matplotlib
- requests


## Setup
Create a local virtual environment and install the project dependencies before opening the code in your editor:

For Windows:
```powershell
py -3.14 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

For Linux and macOS:
```bash
python3.14 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## How to Run
For linux and mac run the following in your terminal:
```
python3 main.py
```
For Windows:
```
py main.py
```
