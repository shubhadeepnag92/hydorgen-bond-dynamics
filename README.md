# hydorgen-bond-dynamics
In this Fortran Program, average number of hydrogen bond, the continuous and the intermittent HB correlation functions are calculated for alpha-crystalline phase of methanol crystal at 90 K.
The correlation functions are calculated from the definitions given in Amalendu Chandra, Phys. Rev. Lett. 85, 768 – Published 24 July 2000
https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.85.768

The fortran program is hbond_num_acf.f95.
A sample trajectory file "npt_500fs.xyz" is given here along with the details of input variables given in "input_HB.dat".

The number of avg. hydrogen is written in "hbnum.dat".
The coorelation functions are written in "hbac.dat" file.

# Program to Compute the Number of Hydrogen Bonds, Continuous HB, and Intermittent HB Correlation Functions

## Overview

This FORTRAN program calculates the number of hydrogen bonds (HB) and computes the continuous and intermittent hydrogen bond correlation functions (C\_{HB}(t) and S\_{HB}(t)) over molecular dynamics (MD) trajectories. The program specifically handles MeOH molecules, computes distances between atoms, checks bonding criteria based on distance and angle, and evaluates HB correlation functions.

## Files

1. **input_HB.dat**: Input file containing the number of trajectories, molecules, atoms, and correlation time steps, as well as the simulation box dimensions.
2. **npt_500fs.xyz**: Input trajectory file containing the coordinates of all molecules at each time step.
3. **npt_sort.xyz**: Output file logging the positions of hydrogen and oxygen atoms for all molecules at each time step.
4. **hblist_OO.dat**: Output file logging molecule pairs where the oxygen-oxygen distance is less than 3.5 Å.
5. **hblist_OH.dat**: Output file logging molecule pairs where the oxygen-hydrogen distance is less than 2.45 Å.
6. **hblist_OOH.dat**: Output file logging molecule pairs where the O-H-O angle is less than 30 degrees.
7. **hblist_mol.dat**: Output file logging the hydrogen bond functions h(t) and H(t) for each molecule at each time step.
8. **hbnum.dat**: Output file recording the average number of hydrogen bonds at each time step.
9. **hbac.dat**: Output file containing the calculated C\_{HB}(t) and S\_{HB}(t) correlation functions.

## Input Variables

- **ntraj**: Number of trajectory steps.
- **nmol**: Number of molecules.
- **natom**: Number of atoms per molecule.
- **ncorr**: Number of correlation time steps.
- **boxlengthx, boxlengthy, boxlengthz**: Dimensions of the simulation box along the x, y, and z directions.

## Program Structure

1. **Initialization**:
   - The input parameters are read from `input_HB.dat`.
   - Variables are allocated dynamically based on the number of atoms and molecules in the system.
   - The output files are opened and initialized with appropriate headers.

2. **Hydrogen Bond Calculation**:
   - For each trajectory frame, the positions of oxygen and hydrogen atoms are read.
   - The distances between oxygen atoms and between hydrogen and oxygen atoms are calculated using the minimum image convention to account for periodic boundary conditions.
   - Hydrogen bonds are identified based on distance criteria (O-O distance ≤ 3.5 Å and O-H distance ≤ 2.45 Å) and angular criteria (O-H-O angle ≤ 30°).
   - The number of hydrogen bonds is counted and logged for each time frame.

3. **Correlation Function Calculation**:
   - Intermittent and continuous hydrogen bond correlation functions (C\_{HB}(t) and S\_{HB}(t)) are computed by checking the persistence of hydrogen bonds over successive time frames.
   - These functions are normalized and written to the output file `hbac.dat`.

4. **Output**:
   - Various output files contain data on atom positions, molecule pairs forming hydrogen bonds, hydrogen bond counts, and the correlation functions.
   - The time taken to run the program is recorded and output at the end.

## Running the Program

1. Prepare the input files (`input_HB.dat` and `npt_500fs.xyz`).
2. Compile and run the program:
   ```
   gfortran hbond_num_acf.f95 -o hbcorrelation
   ./hbcorrelation
   ```
3. The program will output the following files:
   - **npt_sort.xyz**: Contains the positions of oxygen and hydrogen atoms.
   - **hblist_OO.dat**: Logs pairs of molecules where the O-O distance is less than 3.5 Å.
   - **hblist_OH.dat**: Logs pairs of molecules where the O-H distance is less than 2.45 Å.
   - **hblist_OOH.dat**: Logs pairs of molecules where the O-H-O angle is less than 30°.
   - **hblist_mol.dat**: Contains the values of h(t) and H(t) for each molecule.
   - **hbnum.dat**: Records the average number of hydrogen bonds at each time frame.
   - **hbac.dat**: Stores the calculated C\_{HB}(t) and S\_{HB}(t) correlation functions.

## Key Outputs

- **C\_{HB}(t)** and **S\_{HB}(t)**: These functions describe the time-dependent behavior of hydrogen bonds, capturing both intermittent and continuous bond dynamics.
- **Hydrogen bond counts**: The program computes the average number of hydrogen bonds per time frame, which provides insights into the bonding structure of the system.
- **Hydrogen bond data**: Detailed data about which molecule pairs form hydrogen bonds, based on both distance and angle criteria.

## Modifications

To modify the simulation:
- Change the box dimensions or the number of molecules and atoms in `input_HB.dat`.
- Adjust the criteria for hydrogen bond formation by modifying the distance and angle thresholds in the code (default values are 3.5 Å for O-O distance, 2.45 Å for O-H distance, and 30° for O-H-O angle).
- Change the number of correlation time steps by adjusting the value of `ncorr`.

## Contact

For any issues, please contact via email: shubhadeepnag92@gmail.com.
