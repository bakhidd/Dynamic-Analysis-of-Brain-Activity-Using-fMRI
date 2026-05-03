Attention to Visual Motion fMRI data  - GLM and PPI Analyses 
============================================================

This data has been pre-processed using SPM99.
This analysis is in SPM2. Therefore 
ensure 'defaults.analyze.flip' is set to 
1 in the spm_defaults.m file, before proceeding 
further. See http://www.fil.ion.ucl.ac.uk/spm/spm2.html#Compat
for more details

--------------------------------------------------------------------
GLM Analysis:

1. Create a new directory, and copy the conditions.mat file there
2. In matlab (6.0) type
>> spm fmri
3. At the matlab prompt type
>> load conditions
>> whos
4. Press the 'fMRI' button
5. Specify [design]
6. ISI[secs] = 3.22
7. scans per session [360]
8. Specify design in [scans]
9. Select basis set [hrf]
10. model interactions (Volterra) [No]
11. Number of conditions [4]
12. Name for condition 1 [Fixation]
13. vector of onsets - Fixation [fix_index]
13. duration[s] (events=0) [0]
14. Parametric modulation [None]
15. Name for condition 2 [Stationary]
16. vector of onsets - Stationary [stat_index]
17. duration[s] (events=0) [0]
18. Parametric modulation [None]
19. Name for condition 3 [NoAttMot]
20. vector of onsets - NoAttMot [natt_index]
21. duration[s] (events=0) [0]
22. Parametric modulation [None]
23. Name for condition 3 [AttMot]
24. vector of onsets - AttMot [att_index]
25. duration[s] (events=0) [0]
26. Parametric modulation [None]
27. User specified [0]


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
just been updated (you can press the PWD button to return to your
working directory, if necessary).

SPM will now compute the hyperparameters and parameters and 
update the SPM.mat file once again.

--------------------------------------------------------------------------

To interrogate the results press 'RESULTS'.


**** NOTE THAT in the SPMs 'LEFT-IS-RIGHT' **** (see 
earlier comment on spm_defaults.m)

1. Select the SPM.mat file
2. Create a new t-contrast: [0 0 -1 1 0], call it 'Attention'
3. Mask with other contrasts [No]
4. Title for comparison []
5. p value adjustment to control [None]
6. p value [0.001]

You should see Superior Parietal and Dorso-Lateral Prefrontal 
activations, among others.

By selecting overlays->sections, and selecting the normalised structural
image you should be able to identify the anatomy more accurately.

1. Select the SPM.mat file
2. Create a new t-contrast: [-2 0 1 1 0], call it 'Motion'
3. Mask with other contrasts [No]
4. Title for comparison []
5. p value adjustment to control [FWE]
6. p value [0.05]

SELECTING VOIS

1. Go to -9,-99,12
2. Press VOI 
3. Name of region [V2]
4. Adjust data for [effects of interest]
5. VOI definition [sphere]
6. VOI radius(mm) [6]

This saves the relevant info in VOI_V2_1.mat in the working directory

PPI analysis
------------

1. Press PPI, select the previous SPM.mat file
2. Analysis type ? [PsychoPhysiologic interaction]
3. Select VOI_V2_1.mat
4. Include Fixation [No], Stationary [No]
5. Include NoAttMot [Yes], contrast weight [-1]
6. Include AttMot [Yes], contrast weight [1]
7. Name of PPI [V2x(Att-NoAtt)]

This creates the file PPI_V2x(Att-NoAtt).mat. It contains the variable PPI.ppi 
which you can use a user-defined regresser in a new analysis (see below).
See spm_peb_ppi for a description of the PPI data structure.

Create a new directory  and

1. Copy over PPI_V2x(Att-NoAtt).mat
2. In matlab (6.0) type
spm fmri
3. At the matlab prompt type
>> load PPI_V2x(Att-NoAtt)
>> whos
4. Press the 'fMRI' button
5. Specify [design]
6. ISI[secs] = 3.22
7. scans per session [360]
8. Specify design in [scans]
9. Select basis set [hrf]
10. model interactions (Volterra) [No]
11. Number of conditions [0]
12. User specified [3]
13. [360] regressor 1 [PPI.ppi]
14. [360] regressor 2 [PPI.P]
15. [360] regressor 2 [PPI.Y]
16. name of [interaction]
26. name of [Att-NoAtt]
27. name of [V2]

Now

1. Press 'fMRI'
2. Specify [data]
3. Select the newly created SPM.mat
4. Assign the scans from the 'functional' directory
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
2. Create and select a [1 0 0 0] t-contrast
3. Mask with other contrasts [No]
4. Title for comparison [interaction]
5. p value adjustment to control [none]
6. p value [0.01]

The resulting SPM shows areas showing differential connectivity to V2 
due to the attention vs. no attention conditions. The effect in this 
subject is weak.

Will Penny, Lee Harrison and Klaas Stephan, May 2003
