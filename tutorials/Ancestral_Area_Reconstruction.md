# Ancestral Area Reconstruction
## Biogeographic reconstruction with BioGeoBEARs
Elizabeth M. Joyce 

1. Preparing the Ancestral Area Reconstruction
2. Building, running and testing the models with BioGeoBEARS in R
3. Extracting statistics from the Ancestral Area Analysis
4. Conduct a Biogeographic Stochastic Mapping analysis

### Preparing the Ancestral Area Reconstruction

We will run BioGeoBEARS in R, using scripts modified from Nick Matzke's scripts on the [BioGeoBEARs wiki page](http://phylo.wikidot.com/biogeobears). Refer to this page for much more information, and more example scripts and data.

#### Required input
Prior to setting up the model in R, we need to prepare our input files. As a minimum for an ancestral area analysis, we need:
1. Our final chronogram from our [divergence dating analysis](https://github.com/joyceem/MPEP_tutorials/blob/main/tutorials/DivergenceTimeEstimation_BEAST_TargetCapture.md)
	- This is `Meli_3loci_UCLN_fixtree_10M1k-COMBINED-MCC.tre`
2. A geography file `area_codes.txt`. This is a list with all tips in our tree scored for which area they occur in. 

#### Preparing the geographic information
In AAR with the models we are using, we will reconstruct ancestral ranges based on the contemporary range of the extant taxa on the tips of our phylogeny. So, we need to prepare a file where we have scored the extant ranges. 

This geography file needs to be in a specific format for BioGeoBEARS.  Our geography file is called `area_codes.txt`. Open the file in a plain text editor to observe the format. The names of each taxon in our geography file match that in our tree, and their presence and absence from each of the following areas is recorded with 0/1.

For this analysis, we are going to code our taxa into 8 areas:
 ![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdoc5pMAPIbGouHdfyDM7eom0oQLnvHavXorXl09V7nZ4bE5YgbmvNP3NugNKGjT3eo1WYHwBtUFgcmlp544q_EwHwvJxkxCYctjnY7EVgxl0b-rzlX6naTXKgzZffTJNgf2EGzIe2bnCP43vLl-zo0y0Ky?key=Nqw43J4UWJ1wNgJRDTzDDg)

**Score the areas** for all of our tips by determining the distribution of the taxa. You can use books, papers, or online databases such as [POWO](https://powo.science.kew.org/) or [GBIF](https://www.gbif.org/) to do this. Be careful to critically assess the distributions in large online databases though!

We are also going to add an "areas allowed file" `areas_allowed.txt`. This file restricts which combinations of areas are allowed in the analysis. This can restrict the analysis so that it will not come up with area combinations that you might deem unreasonable, and can also help to speed up the analysis. This needs to be thought through carefully, and must make biological sense to your system/group of interest. Open the `areas_allowed.txt` in a text editor to observe its setup.

Which area combinations would you disallow, and why?

#### Stratifying time
Because we are dealing with an ancient biogeographic history of a family over moving continents, we also want to incorporate time stratification into our analysis to change the connectivity of areas over time. There are a number of ways to do this, but we are going to do it by adding:
1. A time periods file time_periods.txt. Open this to see the time periods delimited in our analysis. Why would we choose these periods?
2. A dispersal multipliers file `manual_dispersal_multipliers.txt`. This file is a series of area matrices (1 per time period) with a weighting of how likely dispersal between each area combination is for each time period. Open this in a text editor to observe its format.

### Building, running and testing the models with BioGeoBEARS in R

Load the instance of [R Studio on the workstation](http://10.153.134.10:8787) by going to your internet browser and typing in http://10.153.134.10:8787 (or clicking the link). Enter your login details.

Set the working directory to  `22_ancestral_area_reconstruction`. 
Open the script `1-BGB_DECvDECJ_TimeStr.R` within R Studio. This script will run through the set-up of time-stratified DEC and DEC-J ancestral area analysis. 

This script implements and compares:
1. The standard 2-parameter DEC model implemented in the program LAGRANGE (Ree & Smith 2008); 
2. A DEC+J model implemented in BioGeoBEARS, wherein a third parameter, j, is added, representing the relative per-event weight of founder-event/ jump speciation events at cladogenesis events.  The higher j is, the more probability these events have, and the less probability the standard LAGRANGE cladogenesis events have; 
3. Some standard model-testing (LRT and AIC) is implemented at the end so that users may compare models. 

We will run through the R script `1-BGB_DECvDECJ_TimeStr.R` together line by line. Be sure to pay attention to the input working path and file names.

After you get the DEC model running, you can stop it (it will take longer than we have in the practical to run!). You can load the results of the analysis I have run with this data by entering the following:

```
load("/data_tmp/mpemaster/data/22_ancestral_area_reconstruction/Rdata/MELI_FinalTree_DEC_TimStr.Rdata")

resDEC <- res
```
Now you should have all the parameters you need to build the DEC+J model. Build the DEC+J model following the script. Once you get the DEC+J model running, you can again stop it, and load the results of my previous run:

```
load("/data_tmp/mpemaster/data/22_ancestral_area_reconstruction/Rdata/MELI_FinalTree_DEC+J_TimStr.Rdata")

resDECj <- res
```

You would also normally want to continue to build and test additional models using DIVA-like (Ronquist 1997) and BAYAREA-like models (Landis, Matzke, Moore, & Huelsenbeck, 2013). However, for this tutorial we will only run and compare the DEC and DEC+J models.

At the end of the script, plot the results of the DEC and DEC+J model. 

![](https://github.com/joyceem/MPEP_tutorials/blob/584cb5e001b7e29641fa1b8cd7282336a29f0222/tutorials/images/Pasted%20image%2020250228112904.png)

![](https://github.com/joyceem/MPEP_tutorials/blob/584cb5e001b7e29641fa1b8cd7282336a29f0222/tutorials/images/Pasted%20image%2020250228112956.png)

What do they show? How do they compare?

The advantage of BioGeoBEARs is that it implements many different range evolution models into a likelihood framework, so that we have the option of comparing models to determine the statistically most likely model. 

Look at the model-testing output. You should see something like this:

`restable`
```
	LnL		numparams	d		e		j
DEC	-110.1467574	2		0.01		0.01		0
DEC+J	-100.9219159	3		0.004823255	1.00E-12	0.195466823

```

What do these statistics show? Use the below figure from the [BioGEOBEARS wiki page](http://phylo.wikidot.com/biogeobears) as a hint:
![image](http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_supermodel.png)


`testtable`

![](https://github.com/joyceem/MPEP_tutorials/blob/584cb5e001b7e29641fa1b8cd7282336a29f0222/tutorials/images/Pasted%20image%2020250228124303.png)

Which is the statistically best model - DEC or DEC+J?

Remember that statistical tests also have their caveats and should never be a replacement for common sense and using biologically relevant information. So always interpret any statistical outputs of any software package critically. Also remember that to fit these range evolution models into a likelihood framework, some of the models have to be modified, and differ in their implementation in BioGeoBEARs to how they were originally described. As such, for biogeographic questions always consider which range evolution model and software makes the most sense for your group and data.

## Extracting statistics from BioGeoBEARs analyses

Now we have figures of our ancestral area reconstruction. But if we want to extract exact values and statistics from our analyses, we need to manipulate the model object (e.g. `resDECj`). We will extract the marginal probabilities of each state for each node using script `2-BGB_GetNodeStats.R`.

We will run the script `2-BGB_GetNodeStats.R` together, line by line.

Examine the final table `DECj_range_probabilities.txt`. Along the top row, you have the range, and along the left hand column you have the node numbers. Use the node number labels that we plotted on the tree in `BioGeoBEARS_tree_APEnodelabels.pdf` to interpret this table. 

![](https://github.com/joyceem/MPEP_tutorials/blob/584cb5e001b7e29641fa1b8cd7282336a29f0222/tutorials/images/Pasted%20image%2020250228131249.png)

Extract of `DECj_range_probabilities.txt` table, viewed in Excel:

![](https://github.com/joyceem/MPEP_tutorials/blob/584cb5e001b7e29641fa1b8cd7282336a29f0222/tutorials/images/Pasted%20image%2020250228131444.png)

What is the probability of each range at each node?

## Conduct a Biogeographic Stochastic Mapping analysis

Next, we will conduct some Biogeographic Stochastic Mapping (BSM) using our preferred biogeographic model (i.e. in this case, DEC+J). BSM can be a great complement to a Bayesian/likelihood Ancestral Area Reconstruction, as it can map biogeographic shifts along branches as well as at nodes. Stochastic mapping simulations are simulations conditional on the observed tree and a given model, and the model parameters. If you average many stochastic mapping simulations, you should get the same probabilities of each state that are produced in ancestral states estimation. However, the individual stochastic maps constitute possible histories, and it can be useful to count stochastically mapped events, their dates, etc. For example, it enables us to estimate metrics like the average number of dispersal events per million years between certain areas.

More information about BSM can be found here:
1. Dupin, Julia; Matzke, Nicholas J.; Sarkinen, Tiina; Knapp, Sandra; Olmstead, Richard; Bohs, Lynn; Smith, Stacey (2016).  
[Bayesian estimation of the global biogeographic history of the Solanaceae](http://dx.doi.org/10.1111/jbi.12898). _Journal of Biogeography_, 44(4), 887-899. 
1. Matzke, Nicholas J. (2016). Stochastic mapping under biogeographical models. PhyloWiki BioGeoBEARS website, 2016,  
[http://phylo.wikidot.com/biogeobears#stochastic_mapping](http://phylo.wikidot.com/biogeobears#stochastic_mapping).

Open the R script 3-BGB_BSM.R in R. Execute the script line by line, paying attention to the names of the input files and working paths. This script builds the BSM, runs 100 iterations, and then summarises the 100 iterations, producing a table and histograms with information about the number of anagenetic dispersal, sympatric, vicariance and founder events over time.

Open the first stochastic map, `DECjTimStr_single_stochastic_map_n1.pdf`. 

![](https://github.com/joyceem/MPEP_tutorials/blob/584cb5e001b7e29641fa1b8cd7282336a29f0222/tutorials/images/Pasted%20image%2020250228182437.png)

What does this show?

Explore the summary tables. What does `sumamry_counts_BSMs` show? What does `all_dispersals_counts_fromto_means.txt` show?

By manipulating the BSM output files such as `ana_events_table` further, you can also summarise biogeographic events between areas over time. We won't do this for this tutorial. but as an example, this is the average number of dispersal events between Sunda, Sahul and Wallace over time:

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdjLibLfKDHWqGQVbnQpxoj6JePMwHxAqR9f8liAVKY9H2qAVgQWoordFmGfVzf75t2wyfPIHTd-MvI7uHUtdg4M_7GS5PIiO_QN87Qi3aXd9n3A0jk43dXwqMc4v0_86pV7tAlj3hcHqwW_AWUkJczR_I?key=Nqw43J4UWJ1wNgJRDTzDDg)


For more details about BSM, refer to the BioGeoBEARS wiki page: [http://phylo.wikidot.com/biogeobears#stochastic_mapping](http://phylo.wikidot.com/biogeobears#stochastic_mapping)
