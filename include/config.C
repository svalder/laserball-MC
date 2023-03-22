// ---------------------------------------------------------
// Content: Configuration file used for the laserball diffuser MC
// Author:  Sammy Valder and Martti Nirkko, University of Sussex (2022)
// ---------------------------------------------------------

// Indepedent variable switch: 0 = Light Injection Point, 1 = Glass Bead Concentration
const int VAR = 1;

// Independent variable array
//const double var[] = {7.5, 10.9, 12.1, 15.0, 18.0}; // injection point offsets mm?
const double var[] = {0.5e-3, 1.0e-3, 1.5e-3, 2.0e-3, 2.5e-3}; // density of bubbles in silgel [g/cm³] - TODO use optimum value ~1.57 mg/mL

// Global constants
const int VERBOSE = 1; // verbosity (0-3)
const int NEVENTS = 1e6; // number of photons, standard is 1e6
const int NBINS = 100; // number of bins

// Laser
const double I0 = 1.; // laser intensity, normalised for now
const double NA = 0.2; // numerical aperture for fibre (from R. Ford's thesis)
const double lambda = 4.05e-5; // wavelength [cm]

// Diffuser ball - TODO REVIEW ALL PARAMETERS
const double flaskdiam = 106.; // diameter of the diffuser flask [mm]
const double R = flaskdiam/2.; // radius
const double na = 1.070; // refractive index (air) - TODO is actually 2:1 mixture of SO2:O2 at 1/3 atm!
const double ns = 1.404; // refractive index (SilGel)
const double ng = 1.530; // refractive index (glass)
// SO2 has n = 1.3047 (20 C) from http://www.microkat.gr/msds/Sulfur%20dioxide.html
// SO2 has n = 1.3396 (25 C) from https://pubchem.ncbi.nlm.nih.gov/compound/sulfur_dioxide
// O2 has n = 1.0003 (0 C) from http://www.kayelaby.npl.co.uk/general_physics/2_5/2_5_7.html

// Glass bubbles
const double r_bub = 2.e-3; // median radius of hollow glass microsphere [cm] - 50% percentile
const double d_bub = 8.8e-5; // approx. wall thickness of hollow sphere [cm] - specs from 3M
const double rho_bub = 2.25; // density of borosilicate glass [g/cm³] - specs from 3M
double con_bub = 1.57e-3; // density of bubbles initialised for single use 

// Quartz rod
const double rodlen = 250.; // length [mm]
const double roddiam = 4.; // diameter [mm]
double OFFSET = 0; // light injection point displacement initialised for single use

// Other
const double c0 = 299.792458; // speed of light in vacuum [mm/ns]
const double cs = c0/ns; // speed of light in SilGel [mm/ns]

