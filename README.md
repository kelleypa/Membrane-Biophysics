# Membrane Biophysics
Analyzing GROMACs coarse-grained (CG) molecular dynamics (MD) simulation of model lipid bilayer, determining raft and nonraft domain formations. 

For a more details, check out my Ph.D. Preliminary Exam Report (PreliminaryReport.pdf) in the repository directory.

# Lipid Rafts
Communication between biological cells and with the organelles in a cell is controlled by the membranes that enclose them. They are composed of lipids and proteins. Lipid molecules with a hydrophilic headgroup and a pair of hydrophobic chains form a bilayer roughly 4 nm thick. Proteins are embedded within or are attached peripherally to the bilayer. 

## Lipid-Lipid Interactions

There are many kinds of lipid that differ in molecular structure and affinity for each other. Below illustrates the affinity of saturated, which consists of single bonded chains of carbon, and unsaturated, differentiated by containing double bond(s) =C-C-C= in the lipid chain, and an extremely vital and present lipid: cholesterol.

<p align="center">
  <img width="50%" src="https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/8ebff236-7e70-49c2-a3db-0c0f39fcebee">
</p>

## Domain Formation

These lipid-lipid interactions drives lateral segregation into compositionally distinct domains creating the local environment necessary for the function of resident proteins. 
<p align="center">
  <img width="50%" src="https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/f5f150d8-c8ca-4cfe-9035-dec9d378111d">
</p>

## Lipid Rafts
Lipid rafts are the most studied example. They are domains enriched in sphingolipids and cholesterol molecules tightly packed together that are envisaged as floating in a more loosely packed sea of surrounding lipid. When clustered together, the lipid raft concept posits, the signaling proteins within these nano-sized domains are triggered.

<p align="center">
  <img width="50%" src="https://github.com/kelleypa/Membrane-Biophysics/blob/main/rafts_trimmed_enhanced_reduced.gif">
</p>

# Martini Course-Grained Simulation
We calculate the trajectories of the atoms on the lipid molecules in a bilayer according to Newton’s laws in MD simulations. To allow the study of larger systems over longer timescales, we employed the CG Martini force field [http://cgmartini.nl/].

![martinilipids](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/13d73afd-5c69-4210-a9b1-a16c3c92f686)


# Domain Determination Method
## Probability of Combinations wiht Repetition of 3 different Species of Lipids 
For a fixed number of lipids comprising a leaflet of the membrane, there is a finite number of possible combinations. Thinking statistically as the sliding window is analogously to selecting a set number of colored marbles from a bag, the probability of all possible combinations without repetition is given by:

![image](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/073fc234-1473-4ece-a38d-fd835e18f035)

As an example, let’s take a nonequal mixture of three lipid types denoted as A, B, and C with a number in each leaflet. Either A or B can be thought of a CHOL and SM that have affinity for each other and tend to form domains; C can be thought of a phospholipid. Taking an simple example of only finding three lipids at one time, the possible combinations with corresponding probability is:

![image](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/b79f142c-9c0e-4ef2-9fa8-12c088396484)

The required density threshold of CHOL and SM (A/B) - rich or phospholipid (C)- rich is found by having the hypergeometric probability of all combinations choose r possible lipids inside a sliding window be ≤ 25%. In other words, we selected a 1 in 4 chance of observing a random combination that meets the given threshold for r lipids found in a sliding window as our cut-off value to determine raft-like or nonraft- like domains. Below illustrates the raft-like (red) and nonraft-like (blue) designation.

![image](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/2ea9730d-beb4-4bd6-9f9c-05a953d6ce28)

## Sliding Window
A sliding windows for which the count of individual lipids exceeds the threshold for random mixing with >75% probability were designated PSM-rich/cholesterol-rich and PC- rich depending upon content. As the simulation has periodic boundary conditions, an extension of the boundary of the simulation box allows the sliding window to start at the bottom and slide right; upon reaching the extent of the right, the sliding resets but shifts up one sliding length. Below is an example of a window sliding by a grid length such that each grid box is tallied 9 times.

![DDanimation](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/8121407f-6c4a-41fb-ad1b-078fa39f89bb)

### Voting 
The majority tally of domain type (raft or nonraft) ultimately classifies the region.

![image](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/936ff4e1-3411-478d-9f5f-631a1978878a)

## Window Size 
Lipids were assigned to PSM-rich/cholesterol-rich (raft-like) domains according to the density of their lateral distribution within the plane of each leaflet. The assignment was made by a sliding window method in which a window ~ 2.4 x 2.4 nm (containing on average 11 lipids) in size was tracked ~ 0.8 nm (average separation between lipid molecules) in each direction. A window for which the count of individual lipids exceeds the threshold for random mixing with >75% probability was designated PSM-rich/cholesterol-rich. Windows that do not satisfy this criterion were deemed to be homogeneously mixed (non-raft). The domain type for each ~ 0.8 x 0.8 nm2 area was assigned and tallied 9 separate times by the sliding window, and ultimately classified by majority vote. The lipids inside the area were correspondingly categorized. 

![image](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/b9a04343-fdba-4e7c-92b9-2a1d94ff2053)


# CG Simulation of Effect of Monounsaturated (POPC) vs Polyunsaturated (PDPC) Lipids
### Domains @ 6 μs
Toward understanding the mechanism of action for DHA on lipid raft size relied on a controlled simulation study. This method allowed us to investigate the size of SM-rich/Chol-rich raft-like domains formed in response to DHA and to analyze the associated changes in composition and order of rafts and nonrafts.
The CG simulations began with homogenously mixed bilayers composed of SM/Chol/POPC and SM/Chol/PDPC in 1:1:1 mol ratio, and SM/Chol/POPC/PDPC in 1:1:0.5:0.5 mol ratio. The propensity for SM and Chol to segregate into lipid rafts was then observed over production runs of 6 μs. The figure below shows snapshots of the upper leaflet after 6 μs of simulation for all 3 compositions of membrane. Color-coded circles indicate the lateral location of SM (red), Chol (white), POPC (yellow), and PDPC (blue) molecules, and a color coding of areas indicates the regions identified as raft-like (red) and nonraft (blue), according to the local concentration of lipid. The snapshots illustrate that the formation and size of SM-rich/Chol-rich raft-like domains are enhanced by PDPC. In SM/Chol/POPC, most of the bilayer remains nonraft, and the domains that are SM-rich/Chol-rich are small and few. An increase in the size of SM-rich/Chol-rich domains accompanies the partial replacement of POPC by PDPC in SM/Chol/POPC/PDPC, and the total substitution of POPC with PDPC in SM/Chol/PDPC results in further increase in size.

![image](https://github.com/kelleypa/Membrane-Biophysics/assets/107891103/900b1030-e62d-4b50-b6c7-4acca96e2eb8)

Check out published results: https://www.sciencedirect.com/science/article/pii/S0022316624001743
