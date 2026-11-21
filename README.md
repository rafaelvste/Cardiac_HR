# Cardiac_HR

<p align="center">
  <img src="https://img.shields.io/badge/Language-Fortran-734F96?logo=fortran&logoColor=white" alt="Fortran">
  <img src="https://img.shields.io/badge/Analysis-Nonlinear%20Dynamics-blue" alt="Nonlinear Dynamics">
  <img src="https://img.shields.io/badge/Plots-Gnuplot%20%7C%20Python-orange" alt="Gnuplot and Python">
</p>

The codes in this repository are used to study cardiac behaviors, chaos, periodic structures (i.e., shrimps), and degeneracy in the Hindmarsh-Rose model. They reproduce the results of the following papers:

- **["Biological degeneracy preserves phase diagrams topology in excitable systems"](https://www.researchsquare.com/article/rs-10975481/v1)**, co-authored by Rafael V. Stenzinger and Marcelo H. R. Tragtenberg (2026).

- **["Shrimp hubs in the Hindmarsh-Rose model"](https://doi.org/10.1063/5.0239268)**, co-authored by Rafael V. Stenzinger, Vinícius L. Oliveira and Marcelo H. R. Tragtenberg (2025).

- **["Cardiac behaviors and chaotic arrhythmias in the Hindmarsh-Rose model"](https://doi.org/10.1016/j.chaos.2023.113983)**, co-authored by Rafael V. Stenzinger, Tadeu E. Scalvin, Patrick A. Morelo and Marcelo H. R. Tragtenberg (2023).

**Please cite the respective paper when using any of the codes.**

The Fortran codes can be compiled using GFortran or ifort and include simple Gnuplot scripts to plot the results. Parameters for reproducing the results are provided in the figures, captions, and/or main text of the papers.

---

## Contents

- [Biological degeneracy preserves phase diagrams topology in excitable systems](#1-biological-degeneracy-preserves-phase-diagrams-topology-in-excitable-systems)
- [Shrimp hubs in the Hindmarsh-Rose model](#2-shrimp-hubs-in-the-hindmarsh-rose-model)
- [Cardiac behaviors and chaotic arrhythmias in the Hindmarsh-Rose model](#3-cardiac-behaviors-and-chaotic-arrhythmias-in-the-hindmarsh-rose-model)

---

## 1. Biological degeneracy preserves phase diagrams topology in excitable systems

The manuscript is available on *Research Square* ([https://www.researchsquare.com/article/rs-10975481/v1](https://www.researchsquare.com/article/rs-10975481/v1)).

`KTz_cell_phase_diagram_ISI.f90`, from the previous repository `Reentry_CML_KTz`, reproduces Figs. 1a, 6a, S2a, S2b and S7. For Figs. 1a, 6a and S7, comment out the logistic function and uncomment the hyperbolic tangent function, and vice versa for Figs. S2a and S2b. Replace $x_r$ with $I_e$ (i.e., keep $x_r$ fixed and replace it with $I_e$ in the loop) in this program to obtain Figs. 6a and S2b.

`HR_cell_phase_diagram_ISI.f90` reproduces Figs. 1b, 6b, S6 and S8. Replace $I_e$ with $x_r$ (i.e., keep $I_e$ fixed and replace it with $x_r$ in the loop) in this program to obtain Figs. 6b and S8.

`KTz_cell.f90`, from the previous repository `Reentry_CML_KTz`, reproduces Figs. 2, 7, S3 and S5. Comment out the logistic function and uncomment the hyperbolic tangent function for Figs. 2 and 7, and vice versa for Figs. S3 and S5.

`HR_RK4.f90` reproduces Figs. 3 and 8.

`KTz_cell_bifurc_ISI.f90` reproduces Figs. 4a (top panel) and 9a (top panel).

`KTz_tanh_cell_Lyap_diag_Xr_T.f90` reproduces Figs. 4a (bottom panel) and 9a (bottom panel) when keeping the value of $x_r$ constant.

`HR_cell_bifurc_APD_ISI.f90` reproduces Figs. 4b (top panel) and 9b (top panel).

`HR_cell_Lyap_diag_Ie_b.f90` reproduces Figs. 4b (bottom panel) and 9b (bottom panel) when keeping the value of $I_e$ constant.

Figs. 5a, 5b, S1 and S4 are reproduced by plotting the respective equations.

---

## 2. Shrimp hubs in the Hindmarsh-Rose model

The paper was published in *Chaos* ([https://doi.org/10.1063/5.0239268](https://doi.org/10.1063/5.0239268)).

`HR_cell_bifurc_APD_ISI.f90` calculates the bifurcation diagrams of the ISI in Figs. 1 and 2a.

`HR_cell_Lyap_Ie.f90` is used to calculate the Lyapunov spectrum with respect to the external current parameter $I_e$ in Fig. 2b.

`HR_RK4.f90` is used for the membrane potential results of Fig. 3 and the ECG results in Fig. 9.

`HR_cell_phase_diagram_ISI.f90` is used to generate the HR phase diagrams in Figs. 4, 6, 7 and 8. The Python script `discont_ead_lao.py` determines the discontinuity between maximum ISI values and generates the set of points in Fig. 6b.

`HR_cell_Lyap_diag_Ie_b.f90` is used to calculate the Lyapunov spectrum for different values of $I_e$ and $b$. The results are then processed by the Python script `filter_lyap.py` to combine the first and second exponents into a single file for the Lyapunov phase diagram in Fig. 5. The script `best_fit.py` is used to fit the manually obtained maximum stability points.

The Lyapunov spectrum calculation was adapted directly from the Fortran 77 code in Appendix A of the paper "Determining Lyapunov exponents from a time series" by Wolf et al. ([https://doi.org/10.1016/0167-2789(85)90011-9](https://doi.org/10.1016/0167-2789(85)90011-9)).

### Update — November 2026

A bug in `HR_cell_Lyap_Ie.f90` and `HR_cell_Lyap_diag_Ie_b.f90` affected the results of Figs. 2b and 5. The transient period was not removed when calculating the time average of the sum of the logarithms. The overall appearance of the curve and diagram is correct, but the correct values of the exponents are larger. The bug was fixed in both codes.

---

## 3. Cardiac behaviors and chaotic arrhythmias in the Hindmarsh-Rose model

The paper was published in *Chaos, Solitons & Fractals* ([https://doi.org/10.1016/j.chaos.2023.113983](https://doi.org/10.1016/j.chaos.2023.113983)).

`KTz_cell.f90`, from the previous repository `Reentry_CML_KTz`, reproduces Supplementary Fig. S1.

`HR_cell_bifurc_APD_ISI.f90` calculates the bifurcation diagrams of the action potential duration (APD) and interspike interval (ISI) in Fig. 2. Fig. S2 can be generated by modifying this code with the diastolic interval (DI) condition.

`HR_cell_phase_diagram_ISI.f90` is used to generate the HR phase diagrams in Figs. 3a, 3b, 8a, 8b, 8d, S3, S4, S5, S8 and S9. The hyperbolic KTz phase diagrams in Figs. 3c and 3d can be reproduced using the code `KTz_cell_phase_diagram_ISI.f90` from the `Reentry_CML_KTz` repository by simply replacing the logistic function with the hyperbolic tangent (which is commented out in the code) and setting the parameters according to the figure.

`HR_RK4.f90` is used for the membrane potential results of Figs. 1, 4, 5 and S10. It also allows the cell to be stimulated periodically and the electrocardiogram (ECG) to be generated, which can be used to reproduce Figs. 6, 7, 8c, S6 and S7. The indicated sections of the code responsible for stimulation and ECG generation can be ignored or commented out when generating Figs. 1, 4, 5 and S10.