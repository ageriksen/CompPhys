Readme:
Studies of phase transitions with the Ising model, sans allied magnetic field using Monte Carlo methods with Metropolis Markov chains as sampling rule.

System split into 2 somewhat arbitrary categories, large and small lattice size. Here, large is anything above an axial length of 20.
For the small lattice, variables are controlled against lattice size as well as temperature, regarded as functions of time. Time here, is measured in MC sweeps of the lattice.
For the larger lattice, we express the variables as functions of temperature T around the area of a critical temperature in an attempt to study the phase shift occurring.

Data stored in 2 folders. Large and small L. Large being L>20.
For large L, name arrays as "observable+L##.bin"
For small L, name arrays "observable+L##T##.bin"
T here being temperature and the ## given by:
(int)(tempCount*tempStep + initialTemp)

For plotting and extraction of data,
Files are read in python with glob("filenames") and fed into
Numpy's "fromfile()" function.
Plotting is either matplotlib's plot() or semilogx() as appropriate to task.

MAKE SURE THE FOLDER IN THE FIRST LINE OF THE INPUT FILE EXISTS ALREADY.

=================================
PARALLELLIZATION


Ideas for improvement:
making class objects for
    lattice
        initialize
        hold energy and magnetization
    random
        initializing gen and mt199...
        simple callable
        easy to send into various object in project
    Monte Carlo
        choose random index
        Check with model
        check with sampling rule
        update variables
    Ising model
        model for force interactions between elements.
        Could be in lattice
    Metropolis sampling rule
        Finding rate of probabilities and comparing to variable preference

=================================

