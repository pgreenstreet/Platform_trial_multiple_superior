# Platform_trial_multiple_superior_fixed_allocation_ratio

Here is the code used to produce Table 1 and 2 in the paper, A preplanned multi-stage platform trial for discovering multiple superior treatments with control of FWER and power. 
Here the conjunctivepowerdesign.R is used to create the conjunctive design with the second treatment added halfway through the trial. Similar for Pairwisepowerdesign.R which does this for the pairwise design. 
Disjuctivepowercalculator.R, conjunctivepowercalculator.R and pairwisepowercalculator.R are used to find the powers under each configuration. With conjunctivepowercalculator is only needed when both have a clinically relevant effect. 
Finally expectedsamplesizecalculator.R is used to calculate the expected sample size for each configuration.
