# Long-Term-Cascadia-ETS-velocities
Code that underlies the results of Bartlow (2019).

This code is made avaialble without any promise of user-friendliness.  Cascadia_secular_fits.m is the main code for generating inter-ETS and LTA-ETS velocity fields. The input UNR and PBO datasets are stored in MATLAB structures which may be found on Mendeley Data, link pending.

The simple_inversion_1comp_2019 code performs the inversion for LTA-ETS slip rate on the plate interface. It requires Green's Functions, not computed specifically for this study. The Green's functions used were calculated using PyLith and assuming a 3D model of elastic properties.  See the paper for more details.  

The LTA-ETS velocities and uncertainties used in simple_inversion_1comp_2019 were a combined solution. The combined solution is a weighted average over estimates obtained from both the UNR and PBO datasets with different values for input parameters SSE_L, R_T, and T_min.  This is obtained by running the Cascadia_secular_fits.m code multiple times and combining the results.  More details are given in the paper.
