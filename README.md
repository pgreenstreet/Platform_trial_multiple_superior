# Platform_trial_multiple_superior_fixed_n(k)

Here is the code used to produce the figures in the paper, A preplanned multi-stage platform trial for discovering multiple superior treatments with control of FWER and power. 
Here the conjunctivepowerforfixednk.R is used to create all the sample sizes and boundaries for the conjunctive design. Similar for pairwisepowerforfixednk.R which does this for the pairwise design. One can use this code to also produce the results in table 3 one needs to state how many arms and stages the trial has and when the arms are added. The expected sample size is found using expectedsamplesizefixednk.R.
The two competing single arm designs are given using  singlearmtrialsdesign.R and singlearmtrialsdesignwithFWER.R.
The plots are then produced using plot_maker.R.
