# IncorporationModelMicroMixing
Matlab code and solver for the determination of micro mixing times in a continuous flow setting.
Villermaux-Dushman reaction:
(1) Acid + Buffer <--> Buffer (fast reaction)
(2) Acid + Iodide + Iodate --> Iodine (slow reaction)
(3) Iodine + Iodide <--> Triiodide
This code contains the reaction constants for Perchloric acid (full dissociation), TRIS buffer, and Potassium Iodide/Iodate

## Cite this material
Use DOI:10.5281/zenodo.10103788

## incorporation_model_main.m
Available models: exponential or linear incorporation function (Fournier 1996 or Arian 2021)
Required inputs: Volume flow rates, initial concentrations, experimental measurements of Triiodide concentration & Segregation Index
This file determines the reaction time and micro mixing time (according to Arian 2021)

## ODE_solver_frey.m
The ODE solver returns the triiodide concentration/segregation index for a specific micro mixing time
Model needs to be given

## minimum_I3.m & minimum_xs.m
This function returns the squared difference of the experimental and model Triiodide concentration/segregation index.

## VD_incorporation_model.m
This function returns the temporal evolution of all species concentrations, segregation index and volume flows of the Incorporation Model
