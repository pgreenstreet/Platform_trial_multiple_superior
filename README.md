# Platform_trial_multiple_superior_fixed_allocation_ratio

Here is the code used to produce the Table 1 and 2 in the paper, A preplanned multi-stage platform trial for discovering multiple superior treatments with control of FWER and power. 
Here the conjunctivepowerdesign.R is used to create the conjunctive design with the second treatment added halve way though the trial. Similar for Pairwisepowerdesign.R which does this for the pairwise design. 
Disjuctivepowercalculator.R, conjunctivepowercalculator.R and pairwisepowercalculator.R are used to find the power under each configuation. With conjunctivepowercalculator only needed when both have a clinically relavent effect. 
Finally expectedsamplesizecalculator.R is used to caculate the expected sample size for each configuation.
