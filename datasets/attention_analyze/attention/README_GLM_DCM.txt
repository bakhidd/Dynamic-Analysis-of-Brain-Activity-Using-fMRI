Attention to Visual Motion fMRI data  - GLM-DCM Analyses 
========================================================

This data has been pre-processed using SPM99.
This analysis is in SPM2. Therefore 
ensure 'defaults.analyze.flip' is set to 
1 in the spm_defaults.m file, before proceeding 
further. See http://www.fil.ion.ucl.ac.uk/spm/spm2.html#Compat
for more details

--------------------------------------------------------------------
GLM Analysis:

1. Create a new directory and copy the factors.mat, factor_time_series.mat and corr_look.m files there
2. In matlab (6.0) type
>> spm fmri
3. At the matlab prompt type
>> load factors
>> whos
4. Press the 'fMRI' button
5. Specify [design]
6. ISI[secs] = 3.22
7. scans per session [360]
8. Specify design in [scans]
9. Select basis set [hrf]
10. model interactions (Volterra) [No]
11. Number of conditions [3]
12. Name for condition 1 [Photic]
13. vector of onsets - Photic [phot]
13. duration[s] (events=0) [10]
14. Parametric modulation [None]
15. Name for condition 2 [Motion]
16. vector of onsets - Motion [mot]
17. duration[s] (events=0) [10]
18. Parametric modulation [None]
19. Name for condition 3 [Attention]
20. vector of onsets - Attention [att]
21. duration[s] (events=0) [10]
22. Parametric modulation [None]
23. User specified [0]


SPM now creates an SPM.mat file in your working directory

--------------------------------------------------------------------------

Now assign scans to the design as follows:

1. Press fMRI
2. specify design or data [data]
3. Select the SPM.mat file you've just created.
4. Select 'All' the scans from the 'functional' directory
5. Remove global effects [none]
6. High-pass filter [specify]
7. Cut-off period(secs) [128]
8. Correct for serial correlations [AR(1)]


SPM will now update your SPM.mat file.

We can now estimate the model by pressing the 
'ESTIMATE' button. Then select the SPM.mat file that's
just been updated.

SPM will now compute the hyperparameters and parameters and 
update the SPM.mat file once again.

--------------------------------------------------------------------------

To interrogate the results press 'RESULTS'.


**** NOTE THAT in the SPMs 'LEFT-IS-RIGHT' **** (see 
earlier comment on spm_defaults.m)

1. Select the SPM.mat file
2. Choose the F-contrast for eg. the Motion condition
3. Mask with other contrasts [No]
4. Title for comparison []
5. p value adjustment to control [FWE]
6. p value [0.05]

You should see multiple activations in visual cortex.

By selecting overlays->sections, and selecting the normalised structural
image you should be able to identify the anatomy more accurately.

SELECTING VOIS

1. Go to 0,-93,-3
2. Press VOI 
3. Name of region [V1]
4. Adjust data for [effects of interest]
5. VOI definition [sphere]
6. VOI radius(mm) [6]

This saves the relevant info in VOI_V1_1.mat in the working directory

1. Go to -39,-81,3
2. Press VOI 
3. Name of region [V5]
4. Adjust data for [effects of interest]
5. VOI definition [sphere]
6. VOI radius(mm) [6]

This saves the relevant info in VOI_V5_1.mat in the working directory

Now run corr_look.m to plot the V1 and V5 time series and see how their
correlation changes during attention vs. non-attention blocks

If you now use the F-contrast for eg. the Attention condition you will 
see a number of activated areas.

1. Go to -24,-72,54
2. Press VOI 
3. Name of region [PPC]
4. Adjust data for [effects of interest]
5. VOI definition [sphere]
6. VOI radius(mm) [6]

This saves the relevant info in VOI_PPC_1.mat in the working directory

1. Go to -51,33,24
2. Press VOI 
3. Name of region [PFC]
4. Adjust data for [effects of interest]
5. VOI definition [sphere]
6. VOI radius(mm) [6]

This saves the relevant info in VOI_PFC_1.mat in the working directory

DCM Analysis
------------

These steps will set up the DCM analysis described in

Friston et al. (2003) Dynamic Causal Modelling, Neuroimage, 19(4), 
pages 1273-1302.

1. Press DCM
2. Specify or review [specify]
3. Select the SPM.mat file you created in README_GLM
4. Name for DCM_???.mat [VisFroPar]
5. Select all VOIs in order VOI_V1_1, VOI_V5_1, VOI_PPC_1, VOI_PFC_1
6. Include Photic [Yes]
7. Incluse Motion [Yes]
8. Include Attention [Yes]
9. Make the following intrinsic connections 
   V1 to V5, V5 to V1, V5 to PPC, PPC to V5, PPC to PFC and PFC to PPC
   ie. a hierarchy with bidirectional connections
10. Just connect Photic to V2
11. Connect Motion to modulate the bottom-up connection from V1 to V5
12. Connect Attention to modulate the top-down connections from PFC to PPC and from PFC to PPC.

SPM will now estimate the model.

To look at the estimated model

1. Press DCM
2. Specify or review [review]
3. Select DCM_VisFroPar.mat
4. Threshold [0]
5. You can then eg. display -> Effects of [Photic],[Motion],[Attention] etc.

The numerical results will be slightly different to those in the paper (as 
eg. this is a different subject) but you should still see that 
motion modulates the bottom-up connection between V1 and V5, and 
attention modulates the top-down connections between PPC and V5, and 
PFC and PPC.

Will Penny, Lee Harrison and Klaas Stephan, September 2003
