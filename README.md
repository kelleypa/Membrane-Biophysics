# Membrane Biophysics
Analyzing GROMACs molecular dynamics simulation of model lipid bilayer. 

# Lipid Rafts

Lipid-Lipid Interactions
<p align="center">
  <img width="50%" src="https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/8ebff236-7e70-49c2-a3db-0c0f39fcebee">
</p>


Domain Formation
<p align="center">
  <img width="50%" src="https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/f5f150d8-c8ca-4cfe-9035-dec9d378111d">
</p>

Lipid Rafts

![](https://github.com/kelleypa/Membrane-Biophysics/blob/main/rafts_trimmed_enhanced_reduced.gif)



# Martini Course-Grained Simulation
We calculate the trajectories of the atoms on the lipid molecules in a bilayer according to Newton’s laws in MD simulations. To allow the study of larger systems over longer timescales, we employed the CG Martini force field [http://cgmartini.nl/].

![martinilipids](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/13d73afd-5c69-4210-a9b1-a16c3c92f686)


# Domain Determination Method
## Probability of Combinations wiht Repetition of 3 different Species of Lipids 
We make an assumption that half of all available lipids exist in each leaflet. As an example, let’s take a nonequal mixture of four lipid types denoted as A, B, and C with a number in each leaflet. Either A or B can be thought of a CHOL and SM that have affinity for each other and tend to form domains; C can be thought of a phospholipid. Thinking statistically as the sliding window is analogously to selecting possible colored marbles, the probability of all possible combinations without repetition is given by

![image](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/073fc234-1473-4ece-a38d-fd835e18f035)

![image](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/b79f142c-9c0e-4ef2-9fa8-12c088396484)

The required density threshold of CHOL and SM (A/B) - rich or phospholipid (C)- rich is found by having the hypergeometric probability of all combinations choose r possible lipids inside a sliding window be ≤ 25%. In other words, we selected a 1 in 4 chance of observing a random combination that meets the given threshold for r lipids found in a sliding window as our cut-off value to determine raft-like or nonraft- like domains. Anything that does not meet these two criteria is identified as ’mixed’. Below illustrates the raft-like (red) and nonraft-like (blue) designation.

! [image](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/2ea9730d-beb4-4bd6-9f9c-05a953d6ce28)

## Sliding Window
Windows for which the count of individual lipids exceeds the threshold for random mixing with >75% probability were designated PSM-rich/cholesterol-rich and PC- rich depending upon content.Below is an example of a window sliding by a grid length such that each grid box is tallied 9 times.

![DDanimation](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/8121407f-6c4a-41fb-ad1b-078fa39f89bb)

### Voting 
The majority tally of domain type (raft or nonraft) ultimately classifies the region.

![image](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/936ff4e1-3411-478d-9f5f-631a1978878a)

## Window Size 
Lipids were assigned to PSM-rich/cholesterol-rich (raft-like) domains according to the density of their lateral distribution within the plane of each leaflet. The assignment was made by a sliding window method in which a window ~ 2.4 x 2.4 nm (containing on average 11 lipids) in size was tracked ~ 0.8 nm (average separation between lipid molecules) in each direction. A window for which the count of individual lipids exceeds the threshold for random mixing with >75% probability was designated PSM-rich/cholesterol-rich. Windows that do not satisfy this criterion were deemed to be homogeneously mixed (non-raft). The domain type for each ~ 0.8 x 0.8 nm2 area was assigned and tallied 9 separate times by the sliding window, and ultimately classified by majority vote. The lipids inside the area were correspondingly categorized. 

![image](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/b9a04343-fdba-4e7c-92b9-2a1d94ff2053)


# CG Simulation of Effect of Monounsaturated (POPC) vs Polyunsaturated (PDPC) Lipids
### Domains @ 6 $\mu$s

![image](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/900b1030-e62d-4b50-b6c7-4acca96e2eb8)
