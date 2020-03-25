Patrick's CSB final project2
Description of the project

I am rotating in a Cancer Biology lab that studies the protein Calreticulin (CALR). 
We see that mutant calreticulin drives a cancer phenotype in myeloid cells.
We currently have proteomics data showing all of the proteins mutant calreticulin binds 
in the cell (which is a lot of proteins), but have done nothing with it. The code below, as well
as code run prior (in python) attempts to take those data and overlap them with 15 datasets
identifying genes in different pathways. This way, by the end of this exercise we will have
novel information indicating proteins mutant calreticulin is binding. 
The goals of this project include producing multiple graphs summarizing my findings, as well
as a few summative tables indicating key genes (and their protein products) whose
functions are potentially altered due to mutant calreticulin binding. A general note : most of
the code below was run multiple times--once for each data set. I have tried to group these 
repeated bits of code clearly, so that you can read one line, get the picture, and move on.


FINAL FIGURES
"figure1", "figure2", "figure3", "figure4"
"Table1", "Table2", "Table3"
Compile tables, remove row names##
Table1 = Table shows the most enriched bound proteins (compared to wild type)
Table2 = Table shows the most enriched bound proteins (compared to wild type) which
appear in at least 2 different pathways that i analyzed
Table3 = Table shows the most enriched bound proteins (compared to WT) which appear in at 
least two pathways and excludes all proteasome bound proteins since these are possibly ust
being sent to the proteasome and degraded
figure1 = Figure showing the number of affected genes and their protein products in each pathway
Figure2 = Figure showing the average frequency score of each pathway
Figure3 = Figure showing the percent of each pathway that mutant CALR is binding -- i.e. if
the value on this graph is 100, then 100 percent of proteins in that pathway are being 
bound by mutant CALR
figure4 = dot plot showing the correlation (more accurately the lack of correlation)
between frequency of binding and number of proteins bound in each pathway



