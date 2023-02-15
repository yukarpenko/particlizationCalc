## Instructions to set up the hydrodynamics + polarization toolchain.

Generally, the chain consists of two separate codes: hydrodynamic simulation (vHLLE) which produces the freeze-out/particlization hypersurface,
and this code, which computes the Cooper-Frye-like integrals on the particlization hypersurface in order to compute the polarization of Lambda hyperons or any other hadrons.

### 1. Get the vHLLE code
by following the instructions at:
https://github.com/yukarpenko/vhlle

README.txt, part "I. GETTING and BUILDING vHLLE on Linux"

!! A non-main branch of the code is used to compute polarization. Therefore, after setting up vHLLE, checkout into polarization_dbeta_dec2021 branch. To do that, in vhlle/ directory run `git checkout polarization_dbeta_dec2021` . Re-compile the code with `make clean && make` .

### 2. Getting the 'particlizationCalc' code
get it from this repository:
https://github.com/yukarpenko/particlizationCalc

Clone the repository, cd into the created particlizationCalc/ directory,
checkout into 'polar_xi' branch (`git checkout polar_xi`)

Finally, compile the code:
`mkdir obj; make`  -> which should create a binary named "calc".

**Remarks for Apple users:** \
To compile the code on OSX, some requirements must be satisfied before running `make`. For the following steps it is assumed that Homebrew has already been installed on the system.
1. Natively, the clang compiler does not have access to the OpenMP header file which is needed in the code. To install the library, execute `brew install libomp`. By default, this will create the `omp.h` file in a directory similar to `/opt/homebrew/Cellar/libomp/15.0.7/include/`.
2. To make the OpenMP header file accessible to the compiler, the path to `omp.h` must be set as a CPATH environment variable by running `export CPATH=[path_to_omp.h]`, so for the example path above `export CPATH=/opt/homebrew/Cellar/libomp/15.0.7/include/` (it can be checked that the variable has been set correctly by running `env`)
3. Now, the code can be compiled by running `make` in the particlizationCalc/ directory

### 3. Running the chain
a) run vHLLE for example with: \
`cd vhlle/` \
`./hlle_visc -system RHIC200 -params params/glissRHIC/gliss2RHIC.20-50 -ISinput ic/glissando/sources.RHIC.20-50.dat -outputDir output/rhic200.20-50` \
The hydro simulation above will produce a bunch of output files in the directory specified after the  -outputDir argument. One of the files, which is needed for the next step is "beta.dat"

  b) compute the spectrum+polarization of produced Lambda hyperons \
  `cd particlizationCalc/` \
  `mkdir output` \
  `./calc ../vhlle/output/rhic200.20-50/beta.dat output/rhic200.20-50`
 
### 4. Working with the output
The resulting output file `output/rhic200.20-50` contains a map of numerator and denominator of Eq. 10 in arXiv:1610.04717, in (px,py), or more precisely (pT,phi_p) plane at mid-rapidity. \
The format of the columns is the following: \
col. 0    1     2      3    4    5    6    7     8     9     10 \
     pT phi_p  dN/dpt  s^0  s^1  s^2  s^3  xi^0  xi^1  xi^2  xi^3 \
s^mu is the numerator of Eq. (10) in 1610.04717,  dN/dpt is its denominator - all at a given pT and phi.

In order to compute the transverse momentum-differential polarization, one should essentially divide s^mu by dN/dpt. The repository contains a Python3 script `output/showPolar-xi.py` which does this.
The Python script also does the boost of the polarization vector (in practice, a boost of s^mu) to the rest frame of Lambda -> because it is the quantity which should be compared to experimental data, see Eq. 12 in the same paper.
To run the script type: `cd output; python3 showPolar-xi.py rhic200.20-50`. The argument to the Python script is the relative path to the output file of particlizaitonCalc code. You'll need pyton3, numpy-python3, matplotlib-python3 installed for the script to run.
 
An update: recently, we have published an updated formalism to compute the Lambda polarization: arXiv:2103.14621, which includes an extra term - see Eq. 3 therein.
Therefore, the output from step 3 contains separately the "standard" polarization term from 1610.04717, and also the polarization stemming from the new term - therefore there are columns (xi^0, xi^1, xi^2, xi^3) in the output. So in practice the total polarization would be equal to (s^i + xi^i)/(dN/dpt).
