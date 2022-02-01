# particlizationCalc
This code computes the polarization of (Lambda) particles produced at freeze-out in heavy ion collisions. It works in combination with 
the polarization branches of vhlle.

## Prerequisites
* [vhlle](https://github.com/yukarpenko/vhlle)

* [vhlle_params](https://github.com/yukarpenko/vhlle_params.git)

Once vhlle is correctly installed, be sure to **checkout to one of the polarization branches**. Those that include the "noTgrad" in the name employ the isothermal 
freeze out approximation, for discussion see [2103.14621](https://arxiv.org/abs/2103.14621). 
The polarization branches generate as output a "beta.dat" file, that contains information about the freeze out and is used as input in particlizationCalc.

## Compiling the code

1 clone this repository

2 cd to the particlizationCalc directory

3 create the directory obj: `mkdir obj`

4 compile: `make`

These steps should produce an executable called calc.

## Running the code

The code can be used to compute the polarization of Lambda hyperons at midrapidity. This can be done both for primary and for secondary particles. 
A beta.dat file produced by vhlle is necessary, and in what follows its path will be PATH_TO_BETA. 

### Primary particles
To compute the polarization of primary Lambdas at midrapidity:

`./calc PATH_TO_BETA output_name`

output_name is the name of the output files. This command produces two files:

1 output_name.dim -> contains information about how many different values of transverse momentum pT and azimuthal angle phi_p have been explored.

2 output_name.dat -> a table containing information about the lambda polarization **in the laboratory frame** at different momenta. The table is formatted as:

  pT phi_p  dN/dpt  s^0  s^1  s^2  s^3  xi^0  xi^1  xi^2  xi^3

dN/dpt is the denominator in the mean spin formulae, s^\mu is the numerator of the mean spin induced by thermal vorticity, xi^\mu is the numerator of
the mean spin induced by thermal shear. For the polarization, s and xi must be multiplied by two.

### Lambdas from decays
To compute the polarization transferred from some mother particle to Lambdas at midrapidity: 

`./calc PATH_TO_BETA output_name pdg_mother pdg_second_son`

This command computes the polarization inherited by the Lambda particle in the decay mother -> Lambda + second_son. The arguments pdg_mother and pdg_second_son are 
the Particle Data Group codes associated to the particles (see the [PDG numbering scheme](https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf)).

The decays included up to now are:

* Sigma* -> Lambda + pion (pdg_mother=3224, pdg_second_son=211)
* Sigma0 -> Lambda + photon (pdg_mother=3212, pdg_second_son=22)

The code produces a single file, output_name.dat, containing a table formatted as:

pT phi_p  dN/dpt  tot_S^1  tot_S^2  tot_S^3

where the first three symbols are the same as in the previous paragraph, and tot_S is the total contribution to the numerator of the mean spin vector coming 
by thermal vorticity and thermal shear combined **in the rest frame of the Lambda**. By default, sixteen values of pt (0 to 3 with steps of 0.2) and 
forty values of phi_p (0 to 2π with step π/20) 
are explored.  


