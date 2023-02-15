Files in this folder contain R scripts to reproduce simulation in the manuscript "Theory for Identification and Inference with
Synthetic Controls: A Proximal Causal Inference Framework".

Below are the R scripts to implement the simulation:
* LinearSetting.R: Simulation of the classical setting SC with independent random error across units and time in Section 4 as well
as linear setting with serially correlated error in Section I.1 of the Supplementary Materials  (require myfunctions_LinearSetting.R);
* TimeVaryingEffect.R: Simulation of the linear setting SC with time-varying treatment effects in Section I.2 of the Supplementary Materials 
(require myfunctions_TimeVaryingEffect.R);
* NonlinearSetting.R: Simulation of the Poisson setting in Section I.2 of the Supplementary Materials (require myfunctions_NonlinearSetting.R);
* ConformalLinearSetting.R: Simulation of the conformal prediction intervals based on the proximal inference approach for SC in the linear setting 
in Section A.2 of the Supplementary Materials (require myfunctions_Conformal.R).
