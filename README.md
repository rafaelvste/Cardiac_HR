# Cardiac_HR
The codes in this repository are used in the study of the cardiac behaviors, chaos and periodic structures (i.e. shrimps) in the Hindmarsh-Rose model. They reproduce the results of the papers “Cardiac behaviors and chaotic arrhythmias in the Hindmarsh-Rose model”, co-autored by Rafael V. Stenzinger, Tadeu E. Scalvin, Patrick A. Morelo and Marcelo H. R. Tragtenberg, and "Shrimp hubs in the Hindmarsh-Rose model", co-autored by Rafael V. Stenzinger, Vinícius L. Oliveira and Marcelo H. R. Tragtenberg.

Fortran codes can be compiled using GFortran or ifort and contain simple Gnuplot scripts to plot the results. Parameters to reproduce the results are in the figures, captions and/or in the main text of the papers.

#############################################################################

For “Cardiac behaviors and chaotic arrhythmias in the Hindmarsh-Rose model”:

The paper was published in Chaos, Solitons & Fractals (https://doi.org/10.1016/j.chaos.2023.113983).

KTz_cell.f90, from the previous repository 'Reentry_CML_KTz', reproduces Supplementary Fig. S1.

HR_cell_bifurc_APD_ISI.f90 calculates the bifurcation diagrams of the action potential duration (APD) and interspike interval (ISI) in Fig. 2. Fig. S2 can be generated by modifying this code with the diastolic interval (DI) condition.

HR_cell_phase_diagram_ISI.f90 is used to generate the HR phase diagrams in Figs. 3a, 3b, 8a, 8b, 8d, S3, S4, S5, S8 and S9. The hyperbolic KTz phase diagram of Figs. 3c and 3d can be reproduced using the code KTz_cell_phase_diagram_ISI.f90 from the 'Reentry_CML_KTz' repository by simply changing the logistic function to the hyperbolic tangent (commented in the code) and setting the parameters according to the figure.

HR_RK4.f90 is used for the membrane potential results of Figs. 1, 4, 5 and S10. It also contains the possibility to periodically stimulate the cell and calculate the electrocardiagram (ECG), which can be used to reproduce Figs. 6, 7, 8c, S6 and S7 (the indicated parts of the code responsible for the stimulus and ECG can be ignored or commented out when generating Figs. 1, 4, 5 and S10).

##############################################################################

For "Shrimp hubs in the Hindmarsh-Rose model":

HR_cell_bifurc_APD_ISI.f90 calculates the bifurcation diagrams of the ISI in Figs. 1 and 2a. 

HR_cell_Lyap_Ie.f90 is used for the Lyapunov spectrum with respect to external current parameter Ie in Fig.2b.

HR_RK4.f90 is used for the membrane potential results of Fig. 3 and the ECG results in Fig. 9.

HR_cell_phase_diagram_ISI.f90 is used to generate the HR phase diagrams in Figs. 4, 6, 7 and 8. The Python script discont_ead_lao.py determines the discuntinuity between maximum ISI values and generates the set of points in Fig. 6b.

HR_cell_Lyap_diag_Ie_b.f90 is used for the Lyapunov spectrum for different Ie and b. The results are then processed by the Python script filter_lyap.py to combine the first and second exponents into a single file for the Lyapunov phase diagram in Fig. 5. The script best_fit.py is used to fit the manually obtained maximum stability points.

The Lyapunov spectrum calculation was adapted directly from the Fortran 77 code in Appendix A of the paper "Determining Lyapunov exponents from a time series" by Wolf et al. (https://doi.org/10.1016/0167-2789(85)90011-9).

