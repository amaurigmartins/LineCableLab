[![View OHLToolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/130914-ohltoolbox)

# OHLToolbox

Toolbox in MATLAB for modelling of overhead and underground transmission lines (per-unit-length parameters, propagation characteristics, frequency scan, transient simulations).

### Highlight of the main features

- OHLToolbox is distributed with a GUI designed to aid data entry. You can save your projects in matfiles and recover them for later use.
- Several impedance and admittance formulas, frequency-dependent earth models and modal decompositions techniques are available.
- It is possible to work with overhead, underground and mixed overhead-underground (coming soon) conductor arrangements.
- Parameters are computed for any number of phases with any number of conductor per phase, with or without Kron elimination. 
- You can export the computed data as a JMarti line model (PCH file), compatible with ATP/ATPDraw, or to a EMTP-compliant matfile to be used with the LineCable_Data routine available in the EMTP® software.
- Time-domain transient simulations using the Laplace transform, with several energization sources available.

### Basic instructions

The use of the toolbox is quite self-explanatory. Launch the file 'ohlt.mlapp' from MATLAB workspace and fill in the necessary information. You can save your design by clicking 'Save input session' or recover an existing design by clicking 'Load input session'. 

You will be prompted to specify a JobID and this is **mandatory**. The JobID is a text string which will be used to identify all the output files and folders. When you are satisfied with the data entry, hit the button 'Start'. This will create a subfolder named JobID inside the working directory. Inside this folder you will find two files: 'LineData_fun.m' and 'RunJob_JobID.m'. Call the script 'RunJob_JobID.m' and that's it.

The columns under the header 'OHL cross section & conductor data' are described as follows:

- 1 column -- **PhaseNum**: number of phase (set to 0 for Kron reduction, map to the same phase for bundled conductors).
- 2 column -- **Horiz**: x position of each conduntor in meters.
- 3 column -- **Vert**: y position of each coductor in meters.
- 4 column -- **IntRadius**: internal radii of each conductor in meters.
- 5 column -- **ExtRadius**: external radii of each conductor in meters.
- 6 column -- **CondResist**: resistivity of the conductor in ohms-meter.
- 7 column -- **CondPerm**: relative permeability of the conductor in pu.
- 8 column -- **InsulRadius**: external radii of insulation in meters. Set to NaN if bare conductor.
- 9 column -- **InsulPerm**: relative permeability of insulation in pu. Set to NaN if bare conductor.
- 10 column -- **InsulPermit**: relative permittivity of insulation in pu. Set to NaN if bare conductor.




### Important information

OHLToolbox uses external tools that must be downloaded from the corresponding sources and placed in the appropriate directories, namely:

- VFIT3 - Fast Relaxed Vector Fitting, developed by Gustavsen et al. Available: https://www.sintef.no/globalassets/project/vectfit/vfit3.zip. Unzip the file into the folder 'JMartiModelFun/vfit3'.
- Eigenshuffle - Consistently sorted eigenvalue and eigenvector sequences, developed by John D'Errico. Available: https://www.mathworks.com/matlabcentral/fileexchange/22885-eigenshuffle. Unzip the file into the folder 'mode_decomp_funs'.
- mpmath - A Python library for arbitrary-precision floating-point arithmetic, developed by Fredrik Johansson. Available: https://mpmath.org/.

The Vector Fitting toolbox is used to create the JMarti line models from the modal parameters, if requested by the user. The Eigenshuffle code is used as an additional method to perform the modal decomposition avoiding switchovers in the frequency-domain. The mpmath is used to accuratelly compute Bessel function terms for large arguments.

These codes are properties of the respective authors, with all due credits given. Observe any restrictions and licensing/usage requirements in the websites above.

When calculating the internal impedances of tubular conductors ("pipes"), the native Bessel functions implemented in Matlab may run into overflow issues, yielding +Inf or -Inf results. This is a well-known behavior of besseli(), besselh(), besselj(), besselk(), bessely(), and is due to some floating-point precision shenanigans when the argument of the function is a large number. As a dirty workaround, it is used an external call to a Python script that runs skin effect calculations outside of Matlab. You need to install the mpmath library described above using pip or the method of your preference, and also to inform the full path to your Python interpreter (.exe or binary file) under the tab 'Computation settings' in the GUI.

### Rescrictions of use

We appreciate the interest in our work and we invite the interested users to use our codes as necessary, as long as they are not embedded in any commercial software, which is **strictly prohibited**. However, if you use the OHLToolbox as a part of scientific research, we kindly ask you to refer to our published papers as per the references list on https://www.mathworks.com/matlabcentral/fileexchange/130914-ohltoolbox.
