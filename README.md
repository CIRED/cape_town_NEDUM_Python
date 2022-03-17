# RUN CODE

- Main script is ```main.py```
- To run properly, the code should be put in a folder that also contains the appropriate data, according to the following tree structure (TO BE UPDATED)

```shell
.
├── 2. Data
    ├── 0. Precalculated inputs
    ├── data_Cape_Town
    ├── FATHOM
│   ├── Flood plains - from Claus
│   ├── Land occupation
│   ├── precalculated_transport
│   └── housing_types_grid_sal.xlsx
├── 4. Sorties
└── code capetown python
    ├── _docs
    ├── _old
    ├── _tests
    ├── calibration
    ├── equilibrium
    ├── inputs
    ├── outputs
    ├── LICENSE
    ├── main.py
    ├── manage.py
    ├── plots.py
    ├── README.md
    ├── requirements.txt
    └── setup.py
  
    
```

- In order to integrate packages seamlessly, ```code capetown python``` should be set as a project in Spyder (or any other Python IDE)
- A more extensive user guide is in the process of writing to ensure better reproductability. For now, the code is being reorganized and commented. It will later be optimized and documented as a proper API


# VERSIONING

- The ```main``` branch contains the original code (with some extra features) from the WB working paper (adapted from Matlab to Python): it should run properly
- The ```TLM-edits``` branch contains some code and folder reorganization without rewriting anything: it should run properly
- The ```TLM-write``` branch is the main working branch. It contains some rewriting and commenting: code may not run fully due to ongoing work (progress will be tracked by adding breakpoints to ```main.py``` in Spyder)
- All other branches are for internal use only