import numpy as np

"""
=================================================================================================
                     PROGRAM STRUCTURE

                         ===========
                         | Dataset |
                         ===========
                             |
                          ========
                          | Data |
                          ========
                             |
==============             ======                =============
| Thermostat | ----------- | MD | -------------- | Potential |
==============             ======                =============
                             |
                          ========
                          | main |
                          ========

- The MD class contains the primary logic of the algorithm
- Thermostat and Potential class are separated from it to allow different
thermostats and potentials can be implemented by inheritance and used flexibly.
can be used
- The Data class stores the data stored in the MD simulation and
takes care of the storage
- Dataset is a data set consisting of time, temperature, ...
- Data holds several datasets and some more data that are not time-resolved
stored
- Instead of using getter and setter the members of Data
and Dataset are public, since they are simple data containers.
- main() calls MD with the parameters needed for the task parts
=================================================================================================
"""

# ================================ Potential-class ================================================

# Virtual class from which concrete potentials can be inherited
# (Here only Lennard-Jones necessary, but so you can quickly implement other potentials).
class Potential:
    def V(r2):
        None

    def F(r):
        None

class PotentialLJ(Potential):
    def V(r2):
        # TODO
        # gets squared distance, returns potential
        None
    
    def F(r):
        # TODO
        # gets position returns force
        None

# ------------------------------ End of Potential-class -------------------------------------------

# ================================ Thermostat class ===============================================

# Virtual class from which concrete thermostats can be inherited
class Thermostat:
    def rescale(v, T):
        # gets list of 2d Vectors and temperature
        None

# No thermostat
class NoThermostat(Thermostat):
    def rescale(v, T):
        # gets list of 2d Vectors and temperature
        None

# Isokinetic thermostat for task d)
class IsokinThermostat(Thermostat):
    def rescale(v, T):
        # gets list of 2d Vectors and temperature
        # returns new velocities 
        None

# ------------------------------ End of Thermostat class ------------------------------------------
# ================================ Data-class ===============================================
# One could use 
# `from dataclass import dataclass`
# `@dataclass`
# `class Dataset:`
# Do be more similar to the cpp template.
# Data set for time resolved data
class Dataset:
    def setData(self, t, T, Ekin, Epot, vS):
        self.t = t # double
        self.T = T # double
        self.Ekin = Ekin # double
        self.Epot = Epot # double
        self.vS = vS # 2D Vector

# Return data of the MD simulation
class Data:
    def __init__(self, n, numBins, binSize):
        self.datasets = [] # Time-resolved datasets
        for i in range(n):
            self.datasets.append(Dataset())
        self.rBin = np.zeros(numBins) # Averaged pair correlation function
        self.g = np.zeros(numBins)
        self.r = [] # snapshot of the final position
                    # For task d) it may be useful to use r instead
                    # in the time-resolved datasets instead

    def save(self, fileNameSets, filenameG, filenameR):
        # TODO
        # saves the different data types into files
        None
# ------------------------------ End of Data-Structs ------------------------------------------

class MD:
    def __init__(self, L, N, particlesPerRow, T, potential, thermostat, numBins):
        self.L = L # length
        self.N = N # No particles
        self.potential = potential
        self.thermostat = thermostat
        self.numBins = numBins
        self.Binsize = 0 # TODO

    # Integration without data acquisition for pure equilibration
    def equilibrate(self, dt, n):
        # TODO
        None
        # no return

    def meausure(self, dt, n):
        # TODO
        None
        # returns instance of Data class

    def centerParticles(self):
        # TODO
        None
        # no return

    def calcT(self):
        # TODO
        None
        # returns temperature of current state

    def calcEkin(self):
        # TODO
        None
        # returns kinetic energy

    def calcEpot(self):
        # TODO
        None
        # returns potential energy

    def calcvS(self):
        # TODO
        None
        # returns center of mass velocity 

    def calcDataset(self):
        # TODO
        None
        # returns current Dataset

    def calcDistanceVec(self, i, j):
        # TODO
        None
        # returns distance vector between particles i and j 

    def calcAcc(self):
        # TODO
        None
        # returns accelerations of particles 
# ------------------------------ End of MD-class ------------------------------------------

LJ = PotentialLJ()
noThermo = NoThermostat()
isoThermo = IsokinThermostat()

partPerRow = # TODO
N = # TODO
L = # TODO
numBins = # TODO

# b) Equilibration test
do_b = False
if do_b:
    T = # TODO
    dt = # TODO
    steps = # TODO
    md = MD(L, N, partPerRow, T, LJ, noThermo, numBins)
    md.measure(dt, steps)
    md.save("b_set.txt", "b_g.txt", "b_r.txt")

# c) Pair correlation function
do_c = False
if do_c:
    stringVec = ["0.01", "1", "100"]
    for s in stringVec:
        T = float(s)
        dt = # TODO
        equiSteps = # TODO
        steps = # TODO

        md = MD(L, N, partPerRow, T, LJ, noThermo, numBins)
        md.equilibrate(dt, equiSteps)
        md.meausure(dt, steps)
        md.save(f"c_{T}_set.txt", f"c_{T}_g.txt", f"c_{T}_r.txt")

# d) Thermostat
# TODO