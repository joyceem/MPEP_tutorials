# Divergence time estimation II
## Node dating in BEAST with large datasets
Elizabeth M. Joyce 

1. Introduction
2. Prepare your molecular data for dating
3. Build your BEAST runs in BEAUti
4. Run BEAST
5. Process BEAST output

## Introduction
In this tutorial date our phylogeny of the mahogany family - Meliaceae, using Bayesian inference node dating, implemented in BEAST. Following on from previous days, the supplied data is HybSeq data enriched for the 353 loci of the Angiosperms353 universal bait kit. 

To run this tutorial, you will need to make sure the following programs are downloaded and installed on the workstation:

- [BEAST2](https://www.beast2.org/) (and affiliated programs BEAUti, LogCombiner, and TreeAnnotator)
- [Tracer](https://beast.community/tracer)
- [Sortadate](https://github.com/FePhyFoFum/SortaDate) (and dependencies)
- [AMAS](https://github.com/marekborowiec/AMAS)
- RStudio, with:
	- ape
	- phangorn
	- rBt

It will also be useful to have:
- the BEAST2 suite, 
- Tracer, and 
- FigTree
downloaded onto your local machine.

## Prepare your molecular data for dating
Although Bayesian methods can simultaneously estimate topology and dates of divergence, this can be too computationally intensive for large, NGS datasets and may result in non-convergence (within your lifetime, anyway). There are a number of options for reducing the computational time of a large dataset with Bayesian Inference divergence dating:
	1. Reduce the size of the dataset
		- reduce no. tips
		- reduce no. loci ("gene shopping")
	2. Reduce the complexity of the model
		- e.g. fix tree topology to a robust, previously estimated topology
	
In this tutorial, we will conduct an analysis in BEAST with gene shopping to show how this is achieved. Because we are only using a small subset of the data, we will also fix the tree topology to the robust topology estimated in our ASTRAL analysis.
### Gene shopping
In Bayesian divergence dating methods, the alignment of molecular data is incorporated into the estimation of ages. However, given the computational intensity of Bayesian Inference, it is often not feasible to incorporate hundreds of loci into the analysis. So, the alignment size can be reduced. There are many ways to approach this, but for this exercise, we will use the [SortaDate](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0197433) pipeline.

Before we can work out the "best" genes to use, we must make sure our gene trees are all rooted in the same way, so that we can estimate the tree length, root to tip variance, and bipartition support.

```bash
mkdir -p 20_gene_shopping/00_rerooted_gene_trees

cd 20_gene_shopping/00_rerooted_gene_trees

```
Now, reroot all of the gene trees we previously estimated with IQTREE in `05_MO_fasta_files`. 
To do this, we will reroot using `pxrr`, and we will run this command in a loop so that it reroots all of the gene trees in the `05_MO_fasta_files` directory:
```bash
for tree in ../../05_MO_fasta_files/*.treefile; do pxrr -t "${tree}" -r -g RUTA_Citrus_hystrix,RUTA_Melicope_ternata,RUTA_Ruta_graveolens -o "$(basename "${tree}" .treefile).rr"; done
```

You should see all of your rooted trees in your directory with `ls`. Open a couple in FigTree to check that they are all rooted to our outgroups in the same way.

Now, we can start shopping for our 'best' genes for dating using SortaDate.

Move out of the rerooting directory, into the `20_gene_shopping` directory:
```bash
cd ../
```

Calculate the root-to-tip variance of all the gene trees using the SortaDate script `get_var_length.py`. We can see what information we need for this by running the command:
```bash
python ~/applications/SortaDate/src/get_var_length.py -h
```

Now we can run the command using the following input:
```
python ~/applications/SortaDate/src/get_var_length.py 00_rerooted_gene_trees --flend .rr --outf r2tvar --outg RUTA_Citrus_hystrix
```

and inspect the output:
```bash
ls
head r2tvar
```
Each line represents the statistics for the gene tree, with the first column showing the root-to-tip variance, and the second column showing the tree length.

Now, we can calculate the bipartition support of all gene trees with `get_bp_genetrees.py`:
```bash
python ~/applications/SortaDate/src/get_bp_genetrees.py -h
```
What input do we need?

Now we can run the script with the command:
```
python ~/applications/SortaDate/src/get_bp_genetrees.py 00_rerooted_gene_trees ../06_astral_MO/meliaceae_334_MO_orthologs.ASTRAL.tre --flend .rr --outf bpsupp
```
Inspect the output:
```
ls
head bpsupp
```
What does it show?

Now that we have the root-to-tip variance, tree length and the bipartition support in separate tables (`bpsupp` and `r2tvar`) we need to combine the results together. This can be done with the `combine_results.py` SortaDate script:
```bash
python ~/applications/SortaDate/src/combine_results.py r2tvar bpsupp --outf combined

head combined
```

Now that our statistics are all in one place, we can sort and get the list of the three best genes. We are using three genes for the purposes of this tutorial, but use the right number for your dataset to balance the amount of information included vs. computational time. To choose our 'best' three genes for dating, we will sort the genes according to: 
1. 3 (bipartition support, i.e. similarity to species tree), 
2. 1 (root-to-tip variance), 
3. 2 (tree length). 

```bash
python ~/applications/SortaDate/src/get_good_genes.py -h

python ~/applications/SortaDate/src/get_good_genes.py --max 3 --order 3,1,2 --outf sortadate_3genes_312.txt combined
```

Now look at your results. Which are the best three genes for our dating analysis?
```
head sortadate_3genes_312.txt
```

It should be loci 5913, 6968, and 5333.

Now that we have shopped for the best loci, we can either treat these as separate loci to put into our dating analysis, or concatenate them. For the purposes of this tutorial, we will concatenate the alignments of loci 5913, 6968 and 5333 into one alignment. This can be achieved with the bash script `make_sortadate_alignment.sh.

For this script, you need to include the output file name of the chosen genes (sortadate_3genes_312), the path for the directory with the cleaned alignments, and the file ending for the clean alignments. The script will read the output file of the `get_good_genes.py` script, find the clean alignments in the alignments directory, concatenate the alignments and provide a summary of the concatenated alignment.

You can run the script as follows:
```bash
bash ../script/make_sortadate_alignment.sh sortadate_3genes_312 ../05_MO_fasta_files ortho.aln.clipkit 
```

Inspect the output of the `make_sortadate_alignment.sh` script. Can you find the folder with the individual gene alignments, and the concatenated alignment?

We now have an alignment with three loci `sortadate_3genes_concat_aln.fasta` to proceed with for our dating analysis.

## Build your BEAST runs in BEAUti
Now that we have our reduced alignment to inform the analysis, and our fossil calibrations, we are ready to set up our Bayesian Inference divergence dating analysis. We will put all of this information, plus some additional information (about tree priors, etc), into the BEAUti GUI interface. This GUI interface will help us set up an `.xml` file with all the information for our model, and we will eventually execute the settings in the `.xml` file with BEAST.

On the workstation, up one directory and make a new folder for the dating analysis.

```bash
cd ../
mkdir 21_dating
```


Now, on your **local computer**, open the BEAUti application by clicking on the icon.

### Partitions
Load your alignment in the `Partitions` panel by pressing the `+` button. This would also be the tab to introduce partitions; however, for this tutorial we will analyse all loci together with the same site, clock and tree models. Choose the correct datatype (Nucleotide).

![[Pasted image 20250225144108.png]]
### Site Model
We are not implementing tip dates in this exercise, so skip to the `Site Model` panel next. The Site Model refers to the substitution model applied to your molecular data. We will select a Gamma Site Model for this analysis.
    - Gamma Site Models have 4 gamma categories. This allows for among-site variation in the rate sets.
    - Estimate the Shape at 1.0. This is necessary to divide the Gamma distribution of rates by 4, as set above.
    - Select the substitution model that best fits your data. For this exercise we will select a GTR model. In a GTR model, all pairwise base rates are estimated apart from one (rate CT). Estimate the frequencies in the model.
    
  ![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcy0Npy1Yf87CjUr15vVMmA_hFxt9uZbaBZEPRZR76sC2LjTf9a_5Q5DLpDsnGAiJPzDammi3DeJv5FXW2Fna-kBXYkeHfUgH3hAhHxovMMGgvfMKynM0yugaqWU70cBalHPAQb67e1dcfZBmV4uvd7AUGN?key=Nqw43J4UWJ1wNgJRDTzDDg)

### Clock Model
  Next, in the `Clock Model` panel, we can set up how our molecular clock is applied. Select which type of clock model you would like to use. 
  - We will choose the `Relaxed Clock Log Normal`.
  - Keep the number of discrete rates to -1. This keeps the number of substitution parameters equal to the number of branches in the tree.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfPyUdQVqS6X2d576LnDherqXYiHrKCCGJ_lkYZOjYId-evHgMl3KjcCxrb41VIFH258vaTH8j_3n-yZVZi6GPSefxsCEgGJiX_1UFuqe_YjqnYOIrYdhlVN2roAfoba53DKiTK-knFAYboYbIFzUCTuio?key=Nqw43J4UWJ1wNgJRDTzDDg)

### Priors
In the `Priors` panel we set the priors to feed into our Bayesian Inference model. All priors with a `.s` in the name refer to substitution model priors, those with `.c` refer to clock priors, and `.t` indicate tree priors. The substitution model and clock priors can be changed, but are related to the model settings set up in the `Site Model` and `Clock Model` panels; therefore these should only be changed with good reason. The main priors we need to set in this panel are our tree prior and introduce our fossil calibration priors. 

#### Tree Prior
Select your tree prior, which sets the patterns of speciation and extinction that will be used to model tree generation in the analysis. For this exercise we will choose a Birth Death Model. In practice, it is always a good idea to do a sensitivity analysis to see the effect of your tree models!

You can change the distribution of the `BDBirthRate.t` and `BDDeathRate.t` priors, and their limits to reflect a ‘realistic’ number of speciation and extinction events per year, respectively, if you have this type of information for your lineage. We don't have this sort of information readily available for Meliaceae, and so we will keep these priors at their default values.

![[Pasted image 20250225144552.png]]

#### Fossil calibrations
Yesterday, in [Divergence time estimation I](https://github.com/joyceem/MPEP_tutorials/blob/main/tutorials/DivergenceTimeEstimation_FossilCalibrations.md), we assessed some extinct fossil taxa that have been described for Meliaceae, and their suitability for calibrating nodes. We are going to use the following fossils to calibrate our tree:

| **No.** | **Calibration**                   | **Age (Ma)** | **Node calibrated** | **Reference**                              |
| ------- | --------------------------------- | ------------ | ------------------- | ------------------------------------------ |
| P1      | _†Manchestercarpa vancouverensis_ | 72.1         | crown Meliaceae     | Atkinson (2020)                            |
| P2      | _†Cedrela_ sp.                    | 51           | crown Cedreloideae  | Hickey & Hodges (1975)                     |
| P3      | _†Swietenia miocenica_            | 22.5         | stem _Swietenia_    | Castañeda-Posadas & Cevallos-Ferriz (2007) |

We also want to put a maximum age on the stem of Meliaceae, but we don't have a good fossil for that node. So, we will use a _secondary_ calibration, that is, a calibration taken from a published chronogram, to constrain the analysis:

| **No.** | **Calibration**                                                | **Age (Ma)**          | **Node calibrated** | **Reference**                              |
| ------- | -------------------------------------------------------------- | --------------------- | ------------------- | ------------------------------------------ |
| S1      | RC-complete analysis stem age for stem Meliaceae/stem Rutaceae | 104.68 (96.28-114.03) | stem Meliaceae      | Joyce et al. (2023)                        |

After assessing our fossils, we decided which nodes would be best to calibrate, and based on the topology of our ASTRAL and ML species trees, decided to use the fossils as follows:
![[Pasted image 20250225151413.png]]

Now that we have reliably identified and dated fossils, and know which nodes we want to calibrate with them, we can enter this information into the `Priors` panel. 

Add a new prior for each fossil that will be used for calibration using the above tree and tables as a guide.

To add a prior, click `+ Add Prior` at the bottom of the window.

**Label** the prior with an **informative name without spaces**. Then, **choose the taxa** to be included for this calibration. 
![[Pasted image 20250225151655.png]]
Only check `monophyletic` if you do not want to allow the model to include any other taxa in that clade.

Next, choose the **prior distribution** from the drop-down menu. This is an extremely important decision impacting the results of divergence time estimation and should be considered carefully! For this exercise we will use a uniform distribution with hard boundaries for all of our calibrations. 
![[Pasted image 20250225151830.png]]

Now click the arrow to the left of the prior to set the hard boundaries of the uniform distribution. We will use a secondary calibration for the root of our tree taken from the stem age of Meliaceae/Rutaceae in one of the analyses of Joyce et al. (2023), and use the same HPD interval as estimated in that analysis for the boundaries (96.28-114.03). 
Check `use originate` if you would like to calibrate the stem node, rather than the crown node of your selected taxa. 
![[Pasted image 20250225152040.png]]

For all internal **primary** calibrations, the maximum age of will be set to the maximum age of our root secondary calibration (114.03 Ma), and the minimum age will be set to the age of the relevant fossil.
![[Pasted image 20250225152844.png]]
![[Pasted image 20250225152722.png]]

Repeat these steps to add P2 and P3.
![[Pasted image 20250225153622.png]]
#### Starting tree
Now that our priors are added, we will open the Starting tree panel. Select:
	`View` -> `Show Starting tree panel`

This step is optional, but will allow us to import a tree to use as a starting tree. This will give the tree topology estimation a "head start", as the analysis starts in a topology with a high likelihood and explores the space around this topology rather than wasting time exploring unlikely tree topologies. This will increase the likelihood that the MCMC chains will converge on the globally optimal topology in a reasonable timeframe. 

To give our analysis a robust, likely topology to start with, we want to use the topology estimated in our concatenated analysis with *all* 353 loci. In future, be sure that all the tips in the tree are present in the alignment, and *vice versa*! 

However, BEAST will not run when we use uniform distributions on our calibration priors if the node heights of our starting tree aren’t within the range of the fossil prior distributions. This is because of the hard boundaries of a uniform distribution: if our starting tree and priors don’t fit together at the beginning the analysis won’t be able to start.

There are multiple ways to fix this. One way is to scale the tree using the `Scale` option of the `Starting tree` panel. Look at the branch lengths of your species tree; if scaling by a factor of 10 or 100 increases the branch lengths to fit within the hard boundaries of your calibrated nodes, then it will work. However, if this does not bring your calibrated nodes within range of your prior distributions, a rough, quick dating procedure with a method like Penalised Likelihood can be used. This will make our starting tree ultrametric and to bring the node heights of calibrated nodes to within the range of the prior distributions we set. 

To get your starting tree ready, perform a quick Penalised Likelihood dating procedure in using the `chronos` function of the `ape` package in R. First, open RStudio.

In R, load your packages, set your working directory to your `21_dating folder`, and read your starting tree.
```R
library(ape)  
library(phangorn)  
library(rBt)  

#Be sure to enter the correct path here!
setwd("~/21_dating")

tree <- read.tree("../concat_gene_part.nwk")  
```

Now, visualise the node numbers on your tree, so that you can ascertain the numbers for the nodes that you need to bring within the range of your uniform distribution.

```R
#Visualise node numbers  
plot(tree, cex=1)  
nodelabels(bg=FALSE, frame="n", cex=1, col="red")  
```

![[Pasted image 20250225160753.png]]

Now that you know the node numbers of your nodes with fossil priors (31, 32, 50, 55), we need to make a table listing the node numbers, and the branch length (/age range) that we want those nodes to be within. This can be very narrow, as it is only the starting tree, and when it is uploaded into BEAST the Bayesian analysis will freely estimate the optimal age within the prior distribution set in BEAUti.

```R
#Make calibration file for nodes with hard boundaries 
node <- c(31,32,50,55)  
age.min <- c(103,73,52,24)  
age.max <- c(105,77,55,27)  
soft.bounds <- c(FALSE,FALSE,FALSE,FALSE)  

mycalibration <- data.frame(node, age.min, age.max, soft.bounds) 

#Check table
head(mycalibration)
```

Now run the quick Penalised Likelihood analysis:
```R
tr.dated = chronos(tree, calibration=mycalibration)
```

Don't worry too much if you get the warning message: `false convergence (8)`. Because we are just making a quick starting tree to get our proper BEAST analysis running, the PL analysis does not have to be rigorous - we just need our nodes of the starting tree to be at compatible ages with the priors. 

Check that your nodes with fossil priors are within the range of the uniform prior distribution, and save your starting tree.
```R
#Check branching times on starting tree
branching.times(tr.dated)
ages <- branching.times(tr.dated)
ages <- round(ages,digits=2)

plot.phylo(tr.dated, y.lim=c(-1,length(tr.dated$tip.label)), 
           no.margin = TRUE, cex=1, tip.color=1, edge.color = 1, label.offset=0.5, edge.width=1)
nodelabels(ages,bg="none",frame="none",cex=1,adj=c(1,1),col="tomato",font=2)


#Save tree as newick
write.tree(tr.dated, file="StartingTree.nwk")
```

![[Pasted image 20250225161200.png]]

Have a look at the starting tree. Are the nodes that we have fossils for within the range of the prior distribution we set?

Now we have our starting tree, we need to input it into BEAUti. 

Download your tree from the workstation onto your local computer using PuTTY, CyberDuck or command line.
Now, in a plain text editor, open your starting tree `StartingTree_concat_gene_part.nwk`. Select the tree text (which is in newick format), and copy it. 

In the `Starting tree` panel, select `Newick Tree` from the drop-down menu. Check `Is Labelled Newick`. Then, paste the tree into the `Newick` box.

Keep other settings on default settings.

![[Pasted image 20250225161929.png]]

#### Operators
Because we are only using three genes for the purpose of this tutorial, it is unlikely that BEAST will be able to estimate the 'correct' topology. If you look at the gene trees of loci 5913, 6968 and 5333 generated on previous days, you will see that these genes have a different topology compared to the species tree we generated with all 353 loci. So, for this analysis, we will tell BEAST not to estimate topology based on these three genes (which it usually does at the same time as estimating the ages of the divergences); **we will fix the tree to the same topology as our starting tree**. This also has the advantage of speeding up our analysis, as there are less moves to make at each step. This is not appropriate to do in every case, but there can be sound arguments to fix the tree topology in a Bayesian analysis.

To fix the topology of our tree to be the same as the topology of our starting tree, we will disable the tree rearrangement moves in the `Operators` panel. Select:
	`View` -> `Show Operators panel`

Set the weight of the following operators to 0:
- Subtree Slide (Performs subtree slide rearrangement of tree)
- Narrow Exchange (Narrow exchange performs local rearrangement of tree)
- Wide Exchange (Wide exchange performs global rearrangement of tree)
- Wilson Balding (Performs Wilson-Balding global rearrangement of tree)
![[Pasted image 20250225163807.png]]

Now, the topology of our tree is fixed to our starting tree.

#### MCMC
Finally, we move to the MCMC panel. Here, we will set the parameters for your MCMC chains. 

Set the chain length to 10 million generations.

Keep the “Store Every”, “Pre Burnin” and “Num Initialization Attempts” to their defaults (-1, 0, 10).  

Set the frequency of logging to the `tracelog`, `screenlog` and `treelog` to `Log Every` 1000 generations.

Importantly, change the name of the `tracelog` and `treelog` to something informative, and in a format that enables the tracking of multiple runs, e.g.:
`Meli_3loci_UCLN_fixtree_10M1k-run1.log`
and
`Meli_3loci_UCLN_fixtree_10M1k-run1.trees`

![[Pasted image 20250225164355.png]]

Now we have finished setting up our model in BEAUti! 

Save the model as an `.xml` file with the same name as your `tracelog` and `treelog`, e.g. `Meli_3loci_UCLN_fixtree_10M1k-run1.xml`

![[Pasted image 20250225164607.png]]

## Run BEAST
Now we can run our divergence dating model saved as an `.xml` file in BEAST.

We always want to run multiple MCMC chains of our model to increase our sample size, check convergence and ensure that we have found the global optimum in our tree space. To do this we can edit our `.xml` files directly to quickly set up multiple `.xml` files that can be run simultaneously.
- Open your `Meli_3loci_UCLN_fixtree_10M1k-run1.xml` file in a plain text editor.
- Scroll down to the `tracelog logger` block of the `.xml` file. Change the `fileName` to have a different run name, i.e. `Meli_3loci_UCLN_fixtree_10M1k-run2.log`
![[Pasted image 20250225165235.png]]
- Scroll down to the `treelog logger` block of the `.xml` file. Change the `fileName` to match the new `tracelog` `fileName`, i.e. `Meli_3loci_UCLN_fixtree_10M1k-run2.trees`
![[Pasted image 20250225165340.png]]
- Now save this .xml file under a new name to match the new names of the `treelog` and `tracelog`
	- Save as… `Meli_3loci_UCLN_fixtree_10M1k-run2.xml`

Repeat this with subsequent numbers to create additional BEAST runs. For the purposes of this tutorial, we will run 5 runs, but it's not uncommon to need more.
![[Pasted image 20250225165959.png]]

Having multiple `.xml` files means you can run all `.xml` files at the same time, and the output will be logged with the relevant run name.

Now, transfer your `.xml` files to your `21_dating` folder on the workstation.

On the workstation, navigate to the directory where your `.xml` files are, and execute them using command line. You would normally run them in a way where they cannot be terminated (e.g. with `nohup` or `screen`) because they will take a long time, but for this tutorial we will just run the first few generations of the first `.xml` file. 

```bash
cd 21_dating

~/applications/beast/bin/beast -seed 93 Meli_3loci_UCLN_fixtree_10M1k-run1.xml
```

You should see something like this if it is running:
![[Pasted image 20250225170348.png]]

If you were running all runs at the same time, be sure to select a different starting seed number for each run.

## Process BEAST output

Open the `21_dating` directory from the GitHub Repository. Look at the output from 5 beast runs of our example data. You should see the `.log`, `.trees` and `.state` files from each run.

We can inspect each file by opening them in a plain text editor, or by exploring them visually in Tracer. Tracer also allows us to check the statistics of the runs when we combine the runs, so that we can see whether the runs have arrived at a similar tree ("converged"), and whether we have enough generations to give us some statistical power.

To view the log files visually, open Tracer and load the log files.
- File -> Import trace files
- Select `Meli_3loci_UCLN_fixtree_10M1k-run1`, `2, 3, 4` and  `5.log`
![[Pasted image 20250226163413.png]]
Explore the output of each run: Select different runs, and look at the trace file of all of the different statistics. 

What do you see? What happens if you change the burn in, to discard early generations? What happens to the statistics of the combined file when there is no burn-in? Why?

Adjust the burn-in to 25% (2500000).
- Select the `Combined Trace File` and examine the `Trace log`.
Have the runs converged?

Examine the `Trace Statistics` panel for the `Combined` run. Note that the ESS is >200 for all statistics. This is typically recommended, but this threshold is also somewhat arbitrary, so your results should be considered on a case-by-case basis.

![[Pasted image 20250226164818.png]]

Now that we have explored the convergence and statistics for our runs, let’s combine our runs into one output in LogCombiner.
- Open LogCombiner
- Add the log files from our BEAST runs
- Select the same burn-in that we decided was appropriate from examining the Tracer traces (25%)
- Check `Resample states at a lower frequency`. We’ll resample every 10000 states to reduce the size of the combined log file.
- Select a location for the `Output File`, and name it something informative (e.g. `Meli_3loci_UCLN_starttree_10M1k-COMBINED.log`.

![[Pasted image 20250226165243.png]]

 Select `Run`.
![[Pasted image 20250226165416.png]]

Repeat the above to combine the tree files.

Open the combined `.log` file and `.trees` file to inspect the output.

Now we have our runs subsampled and combined into one log and tree file, we need to summarise the information from all of these trees together into one, final tree with the posterior probabilities of the nodes, the posterior estimates and HPD limits of the node heights and (in the case of a relaxed molecular clock model) the rates. We will do this in TreeAnnotator.

Open TreeAnnotator
- Select your combined tree file as the `Input Tree File`
- We don’t want to discard any trees in the burn-in because we already did this when we combined the runs: leave the `Burnin percentage` at 0.
- Leave the `Posterior probability` limit to 0.0, as we want PPs to be estimated for all nodes.
- Select `Maximum clade credibility tree` as the `Target tree type`.
	- This will summarise information on the tree in the posterior sample that has the maximum sum of posterior probabilities on its n − 2 internal nodes. This tree is not necessarily the majority-rule consensus tree. If you select the “User target tree” then the tree statistics will be summarised on a user-specified tree.
- Select `Mean heights` for the `Node heights` option.
	- Here you can choose how node heights (divergence ages) are summarised. You can choose to keep the heights that the target tree has (`Keep target heights`), or rescale it to reflect the posterior mean/median node heights for clades.
- Select a location for the `Output File`, and name it something informative (e.g. `Meli_3loci_UCLN_fixtree_10M1k-COMBINED-MCC.tree`.
- Check the `Low memory` box.
![[Pasted image 20250226170659.png]]
Select Run

![[Pasted image 20250226170709.png]]

Now, we have our final chronogram with divergence ages for all our taxa!

Open your chronogram in FigTree to see your final, dated tree (or import and plot in R)!

What does it show? Can you see how it has changed from the starting tree?
![[Pasted image 20250226171248.png]]
