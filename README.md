# Probabilistic Time Petri Nets for Clinical Pathway Design and Analysis
A performance bound solver for Probabilistic Time Petri Nets

The solver computes:
1. Transitions maximum throughput bound 
2. Transitions minimum cycle times (in case of live nets) and identifies the slowest subnet.

The computation is performed using Linear Programming techniques.

## Structure of the repository

- *doc*: documentation (see [Solver modules](https://github.com/simber72/PTPNperfbound/blob/main/doc/Solver_modules.md))
- *examples*: PTPN model examples and results
- *PTPN*: source files for PTPN
- *src*: source files (solver modules)
- *tests*: source files (test scripts)
- *LICENSE*: license file
- *README.md*: this file
- *README.txt*: short explanation of the solver
- *requirements.txt*: Python libraries required to run the solver

## Dependences
- ```cplex```: the optimization problems are solved using IBM CPLEX Studio or CPLEX Studio Community (up to 1000 variables/ contraints). 

Besides other Python external libraries are needed, which can be installed with
 ```pip```: all the dependencies are included in the  ```requirements.txt```).

## How to install

The solver can currently run with Python3 (versions 3.8/3.9/3.10).

0. Clone or download this repository, and move into.

1. Install a virtual environment of Python3 (one of the versions indicated above):

 ```$ python3 -m venv venv``` ...and activate it: ```$ source venv/bin/activate```

2. Update the pip and setuptools versions:

 ```python3 -m pip install --upgrade pip setuptools```

3. Install the libraries required to run the solver:

 ```$ pip install -r requirements.txt``` 

4. Set the Python path environment variable  to the current path:

 ```export PYTHONPATH=.```

## How to use
The soler can be run using a Command Line Inferface (CLI):

```
python3 src/ptpnbound.py --help
Usage: ptpnbound [OPTIONS] NAME TNAME

Options:
  -lp, --lpmodel TEXT          LP model files (CPLEX  models - lp format)
  -lpo, --lpoutput TEXT        Result files (CPLEX model results - xml format)
  -o, --output <TEXT TEXT>     Result file: <name format> (available formats:
                               pnml, dot)
  -v, --verbose                Print results to stdin
  --help                       Show this message and exit.
```
where ```NAME``` is the pathname of the PTPN model (.pnml) and ```TNAME``` is the name of the 
transition of reference for the bound computation. 

In case the last installation step (5.) has not been performed, you can launch the CLI as follows:

 ```python3 src/ptpnbound.py --help```

## Authors
- Manon Le Moigne (ENS Paris-Saclay)
- Cristian Mahulea (Aragon Institute for Engineering Research (I3A), University of Zaragoza)
- Grégory Faraut (Lurpa, ENS Paris-Saclay)
- Simona Bernardi (Aragon Institute for Engineering Research (I3A), University of Zaragoza)
- Jorge Albareda (Hospital “Lozano Blesa” in Zaragoza)
- Lidia Castán (Hospital “Lozano Blesa” in Zaragoza)

## License
GNU GPLv3

## References
M. Le Moigne, C. Mahulea, G. Faraut, S. Bernardi, J. Albareda, L. Castan, "Probabilistic Timed Petri Nets for Clinical Pathway Design and Analysis: A Case Study," Discrete Event Dynamic Systems, vol. 35, pp. 205–231, September 2025. DOI: 10.1007/s10626-025-00419-4.

S. Bernardi, J. Campos, "A min-max problem for the computation of the cycle time lower bound in interval-based Time Petri Nets," IEEE Transactions on Systems, Man, and Cybernetics: Systems, 43(5), September 2013.

Y. Emzivat, B. Delahaye, D. Lime, O.H. Roux, "Probabilistic Time Petri Nets", HAL Open Science,
https://hal.science/hal-01590900.

M. Le Moigne, "ARPE clinical pathway identification algorithms", GitLab repository: https://gitlab.crans.org/bleizi/arpe-identification

S. Bernardi, "PTPNperfbound" ,Github repository: https://github.com/simber72/PTPNperfbound

WSU CASAS smart home project: D. Cook, A. Crandall, B. Thomas, and N. Krishnan. CASAS: A smart home in a box. IEEE Computer, 2013. http://eecs.wsu.edu/~cook/pubs/computer12.pdf

## Funding
This work was supported by the Spanish Ministry of Science and Innovation through the project TED2021-130449B-I00 and by the Aragonese Government under Programa de Proyectos Estratégicos de Grupos de Investigación (COSMOS research group, ref. T64-23R).

<p align="center">
  <img src="./docs/logo_proy.jpg" alt="Logo del proyecto" width="200"/>
</p>
