[![View OHLToolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/130914-ohltoolbox)

# OHLToolbox

Toolbox in MATLAB for modelling of overhead and underground transmission lines (per-unit-length parameters, propagation characteristics, frequency scan, transient simulations).

### Highlight of the main features

- OHLToolbox is distributed with a GUI designed to aid data entry. You can save your projects in matfiles and recover them for later use.
- Several impedance and admittance formulas, frequency-dependent soil models and modal decomposition techniques are available.
- It is possible to work with overhead, underground and mixed overhead-underground (coming soon) conductor arrangements.
- Parameters are computed for any number of phases with any number of conductors per phase, with or without Kron elimination. 
- You can export the computed data as a JMarti line model (PCH file), compatible with ATP/ATPDraw, or to a EMTP-compliant matfile to be used with the LineCable_Data routine available in the EMTP® software.
- Time-domain transient simulations using the Laplace transform, with several energization sources available.

### Basic instructions

The use of the toolbox is quite self-explanatory. Launch the file 'ohlt.mlapp' from MATLAB workspace and input the necessary information. You can save your design by clicking 'Save input session' or recover an existing design by clicking 'Load input session'. 

You will be prompted to specify a JobID and this is **mandatory**. The JobID is a text string which will be used to identify all the output files and folders. When you are satisfied with the data entry, hit the button 'Start'. This will create a subfolder named JobID inside the working directory. Inside this folder you will find two files: 'LineData_fun.m' and 'RunJob_JobID.m'. Call the script 'RunJob_JobID.m' and that's it. Note that the GUI does not actually perform any calculations, it is only a wrapper to the main function call.

The columns under the header 'OHL cross section & conductor data' are described as follows:

- 1 column -- **PhaseNum**: number of phase (set to 0 for Kron reduction, map to the same phase for bundled conductors).
- 2 column -- **Horiz**: y position of each conductor in meters.
- 3 column -- **Vert**: z position of each conductor in meters.
- 4 column -- **IntRadius**: internal radius of each conductor in meters.
- 5 column -- **ExtRadius**: external radius of each conductor in meters.
- 6 column -- **CondResist**: resistivity of the conductor in ohms-meter.
- 7 column -- **CondPerm**: relative permeability of the conductor in pu.
- 8 column -- **InsulRadius**: external radius of insulation in meters. Set to NaN if bare conductor.
- 9 column -- **InsulPerm**: relative permeability of insulation in pu. Set to NaN if bare conductor.
- 10 column -- **InsulPermit**: relative permittivity of insulation in pu. Set to NaN if bare conductor.

### Acknowledgements

The OHLToolbox has been actively maintained and improved by Theofilos Papadopoulos (thpapa@gmail.com), Andreas Chrysochos (anchryso@gmail.com) and Amauri Martins-Britto (amaurigmartins@gmail.com). The authors acknowledge and wish to thank for the efforts of the external contributors, namely: Caio Moraes (implemented the original bundle reduction codes) and Jaimis Leon (provided invaluable insights that enabled us to troubleshoot the VF solution when building the PCH file).

### Important information

OHLToolbox uses external tools that must be downloaded from the corresponding sources and placed in the appropriate directories, listed below:

- vfit3 - Fast Relaxed Vector Fitting, developed by Gustavsen et al. Available: https://www.sintef.no/globalassets/project/vectfit/vfit3.zip. Unzip the file into the folder 'JMartiModelFun/vfit3'.
- Eigenshuffle - Consistently sorted eigenvalue and eigenvector sequences, developed by John D'Errico. Available: https://www.mathworks.com/matlabcentral/fileexchange/22885-eigenshuffle. Unzip the file into the folder 'mode_decomp_funs'.
- mpmath - A Python library for arbitrary-precision floating-point arithmetic, developed by Fredrik Johansson. Available: https://mpmath.org/. Installation details below.
- Bodefit - Bode process program, developed by E. S. Bañuelos-Cabral et al. Available: https://www.intechopen.com/chapters/57667. Deploy the file 'Bode_process.m' to the folder 'JMartiModelFun/functions'.
- Carson's integral - Closed-form solution proposed by T. P. Theodoulidis. Available: https://www.mathworks.com/matlabcentral/fileexchange/50134-carson-s-integral. The files are packaged in the toolbox and no further action is required. 
  

The vector fitting toolbox is used to create the JMarti line models from the modal parameters, if requested by the user. In case the VF does not converge to a stable solution (real poles only), the Bodefit is used as alternative, at the cost of accuracy. The Eigenshuffle code is used as an additional method to perform the modal decomposition avoiding switchovers in the frequency-domain. The mpmath library is used to accuratelly compute Bessel function terms for large arguments.  
These codes are properties of the respective authors, with all due credits given. Observe any restrictions and licensing/usage requirements in the websites above.

Note about Bessel corrections: when calculating the internal impedances of tubular conductors ("pipes"), the native Bessel functions implemented in Matlab may run into overflow issues, yielding +Inf or -Inf results. This is a well-known behavior of besseli() and besselk(), and is due to some floating-point precision shenanigans when the argument of the function is a large number. As a dirty workaround, it is used an external call to a Python script that runs skin effect calculations outside of MATLAB. You need to install the mpmath library described above using pip or the method of your preference, and also to inform the full path to your Python interpreter (.exe or binary file) under the tab 'Computation settings' in the GUI. This is not mandatory, but it may come necessary if you work with hollow tubular conductors.

### On the computational methods implemented in OHLToolbox

The toolbox incorporates several impedance and admittance formulas, under different assumptions to represent the imperfect earth in the line models. For a compreheensive analysis, please refer to the following scientific publications:

- [Closed-form expressions for the analysis of wave propagation in overhead distribution lines](https://www.mdpi.com/1996-1073/13/17/4519)
- [A generalized model for the calculation of the impedances and admittances of overhead power lines above stratified earth](https://www.sciencedirect.com/science/article/pii/S0378779610000684?via%3Dihub)


### Restrictions of use

We appreciate the interest in our work and we invite the interested users to use our codes as necessary, as long as they are not embedded in any commercial software, which is **strictly prohibited**. However, if you use the OHLToolbox as a part of scientific research, we kindly ask you to refer to our published papers:

- T.A. Papadopoulos, Z.G. Datsios, A.I. Chrysochos, A.G. Martins-Britto, G.K. Papagiannis, "Transient induced voltages on aboveground pipelines parallel to overhead transmission lines," Electric Power Systems Research, 109631, ISSN 0378-7796, Jun 2023, doi: 10.1016/j.epsr.2023.109631.

- T. A. Papadopoulos, A. I. Chrysochos, C. K. Traianos, and G. Papagiannis, “Closed-Form Expressions for the Analysis of Wave Propagation in Overhead Distribution Lines,” Energies, vol. 13, no. 17, p. 4519, Sep. 2020, doi: 10.3390/en13174519.

- A. I. Chrysochos, T. A. Papadopoulos and G. K. Papagiannis, "Robust Calculation of Frequency-Dependent Transmission-Line Transformation Matrices Using the Levenberg–Marquardt Method," in IEEE Transactions on Power Delivery, vol. 29, no. 4, pp. 1621-1629, Aug. 2014, doi: 10.1109/TPWRD.2013.2284504.

- A. G. Martins-Britto, T. A. Papadopoulos, Z. G. Datsios, A. I. Chrysochos and G. K. Papagiannis, "Influence of Lossy Ground on High-Frequency Induced Voltages on Aboveground Pipelines by Nearby Overhead Transmission Lines," in IEEE Transactions on Electromagnetic Compatibility, vol. 64, no. 6, pp. 2273-2282, Dec. 2022, doi: 10.1109/TEMC.2022.3201874.
