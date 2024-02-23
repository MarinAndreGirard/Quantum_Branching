###
Jupyter notebooks (.ipynb) worth looking at:
-Basic_tools
-Neff_study
-overlap_study
-distrib_video (or simply look at its product, the pre-generated gifs)

In progress/exploratory notebooks:
-Energy non-conservation in branching
-Weight_conservation

The python code (.py) contains functions used repeatedly.

Put explanation in each notebook.

--Basic tools: Has some graphs to verify that a measurement is indeed happening

--Neff_study: Looks at the distribution of the global state over the total hamiltonian spectrum as well as how Neff_schmidt/Neff_total evolves.
We look at various interaction strengths, which vary the distribution of teh global state over the energy eigenstates. 
In all cases with interaction (and w=0.5) the curves converge to the same value and oscilate around it.
As EI increases, the initial value of Neff_schmidt2 diminishs, with the extrem case of EI = 0.2 making Neff_schmidt 1 and 2 equal and EI=0.25 making Neff_schmidt 2 < Neff_schmidt 1 at the starting point. As EI increases the value at which teh curves converge seems to be diminishing.
We also have a a graph attempting to reproduce one of Andreas graphs, looking at a longer time scale for EI = 0.09
We also look at what varying w between 0.5 and 0.1 does. We see that w=0.5 has the lowest Neff for the gobal state with w=0.3 the largest and it being still higher for 0.2 and 0.1.
Changing w also creats interstig behavior in Neff_schmidt/Neff_tot. Weirdly for w= 0.5 the curves no longer converge. The curves are essentially matched for w=0.3. And Neff_schmidt2/Neff_tot seems less turbulent than schmidt 1 for w = 0.1,0.2,0.3.
There is also a plot of Neff_schmidt/Neff_tot for w=0.3 and t=100 (long time scale) We see that the 2 curves do stabilze to different values.

--overlap_study: Look at the overlap of Schmidt states with each other as well as with the gobal state. Since the schmidt are orthogonal we look at projections of the probability vector (sqrt of it), the inner product wouldd not work.
The first set of graphs look at the impact of w, for 0.5, 0.4, 0.3. We see that the global-schmidt overlap for schmidt 1 and 2 converge closer for w = 0.5 than w=0.3. We also see that for w=0.5 the schmidt 1 and 2 overlap is very noisy for w=0.3,0.4 but not w = 0.5. For w = 0.4, 0.3, the schmidt 1 and 2 overlap starts at the same value as global and schmidt 2 as expected, then dips low, before going back up to a value lower than teh other 2 overlaps and oscilating there.




