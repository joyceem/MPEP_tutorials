# Ancestral State Reconstruction in R
Elizabeth M. Joyce

## Introduction
As discussed during the lecture, there are many different ways to model the evolution of traits. For out practical, we are going to use a simple model to reconstruct the ancestral states of two traits: seed type, and sexual system. Note that these are **categorical traits**, i.e., all of the states are categories (not continuous). There are also multiple options for modelling the evolution of continuous traits, but we will focus on categorical traits for this tutorial.

## Scoring traits

### Seed type in Meliaceae
There are two major seed types in Meliaceae: **winged seeds**, whereby the testa is thin and dry and forms a wing at one end or the whole way around the seed; and seeds without wings (herein referred to as '**unwinged seeds**') whereby the testa is unwinged and either forms a sarcotesta or is covered with an aril.

![](https://github.com/joyceem/MPEP_tutorials/blob/33c6785e67250bf621668255f5ad11c4eb4569af/tutorials/images/Pasted%20image%2020250309105859.png)

A) An example of **winged seeds** in *Toona sinensis* (image: [Roger Culos](https://de.wikipedia.org/wiki/Toona_sinensis#/media/Datei:Toona_sinensis_MHNT.BOT.2010.12.11.jpg)); B) An example of **unwinged seeds** in *Aglaia meridionalis* (image: Elizabeth Joyce), whereby the testa is covered by a bright orange aril.

Did winged seeds evolve once or multiple times? What was the ancestral state of the family? Did winged seeds evolve from unwinged seeds, or vice versa? How frequently do transitions between these states occur? Is it more frequent for winged seeds to evolve into unwinged seeds? Is one state linked to more diversification than the other? These are the sorts of questions that we can tackle with different types of ASR models.

To conduct an ASR, we first need to gather all the information about what we know about seed states in the family. A simple way to do this is by scoring the states of the seed trait for all extant lineages we have in our phylogeny (our tips).

Seeing as there are two potential states, this is a **binary trait**. We will code states of our tips using the following system:
 - 0=unwinged
 - 1= winged
	
Download the `meli_traits_template.xlsx` and use the last generic monograph of the family by Pennington and Styles (1975) to score the states for each of the tips.
### Sexual systems in Meliaceae
There is a great amount of variation in the sexual systems of species in Meliaceae. Members of Meliaceae can be have unisexual flowers in the same individual (monoecious plants); they can have unisexual flowers in distinct individuals (dioecious plants); they can have unisexual and bisexual flowers (polygamous plants), or they can have only bisexual flowers (hermaphroditic plants). This can vary between species within a genus, and sometimes possibly even within a species (but needs further research).

![](https://github.com/joyceem/MPEP_tutorials/blob/90a90e886ffa5f455e03c8b6015a51c2c8228c94/tutorials/images/Pasted%20image%2020250309113552.png)
Modified from https://commons.wikimedia.org/wiki/File:Monoecy_dioecy_en.svg.

How did these different sexual systems evolve?

Seeing as there are four potential states, this is a **multistate trait**. We will code states of our tips using the following system:
 - 0=unisexual flowers in the same individual (monoecious plants)
 - 1=unisexual flowers in distinct individuals (dioecious plants)
 - 2=unisexual and bisexual flowers (polygamous plants)
 - 3=only bisexual flowers (hermaphroditic plants)
	
Download the `meli_traits_template.xlsx` and use the paper: Laino Gama R, Muellner-Riehl AN, Demarco D, Pirani JR (2021) Evolution of reproductive traits in the mahagony family (Meliaceae). _Journal of Systematics and Evolution_ **59**(1), 21–43. [https://doi.org/10.1111/jse.12572](https://doi.org/10.1111/jse.12572) to score the states for each of the tips.

### Formatting your trait matrix
Different ASR methods will want different formats for your trait matrix. Make sure you understand what the requirements are for the ASR you want to run, and format your table accordingly. I have already done this for you

## Modelling ASR in R
There are different programs you can use to run different types of models (e.g. RevBayes, R, Mesquite, etc). We are going to use R for this tutorial.

Load the instance of [R Studio on the workstation](http://10.153.134.10:8787) by going to your internet browser and typing in http://10.153.134.10:8787 (or clicking the link). Enter your login details.

Set the working directory to  `23_ancestral_trait_reconstruction`. Open the script `ASR_MPEP.R` within R Studio. 

### Evolution of discrete characters with Markov Models
To begin, we will model trait evolution under a simple Markov Model, where the rate class is the same across the whole tree. To do this, we will use the R package `corHMM`. For more information on `corHMM` see the manual for the package: https://cran.r-project.org/web/packages/corHMM/corHMM.pdf, and the paper that describes it: Beaulieu et al. (2013) *Systematic Biology* [doi.org/10.1093/sysbio/syt034](https://doi.org/10.1093/sysbio/syt034).

Read in the chronogram we generated earlier in the week, and the trait data matrix. Check that everything is read in properly, and that all of your taxa in your trait matrix are represented in the tips of your tree.

```
## Read in your tree file

tree <- read.tree("data/Meli_chronogram.nwk")

## Check your tree and the format

tree #tree properties
tree$tip.label 
tree$edge #this is how the tree is stored; relationship between nodes and their descendants
tree$edge.length #this is the length of the branch lengths
plot(tree)

## Now load your trait matrix with states scored for each trait

data <- read.csv("data/meli_traits.csv", header=F) # Matrix of seed type and sexual system

## Check the format of the table

head(data)
  # V1=species
  # V2=seed type: 0=unwinged/1=winged
  # V3=sexual system: 0=unisexual flowers in the same individual (monoecious plants); 1=unisexual flowers in distinct individuals (dioecious plants); 2=unisexual and bisexual flowers (polygamous plants); 3= only bisexual flowers (hermaphroditic plants)
  # polymorphisms handled with "&"
tail(data)
```
#### Modelling seed evolution (binary trait)
Now, we are going to compare the results of an ER (Equal Rates) vs. ARD (All Rates Different) model of evolution for our binary seed trait.

![](https://github.com/joyceem/MPEP_tutorials/blob/90a90e886ffa5f455e03c8b6015a51c2c8228c94/tutorials/images/Pasted%20image%2020250309120446.png)

Run your ER model:
```
ans <- rayDISC(tree, data, ntraits=1, charnum=1, model="ER")
```

Look at the results of the model. 
```
> ans
Fit
      -lnL      AIC     AICc ntax
 -11.05538 24.11076 24.27076   27

Rates
            0           1
0          NA 0.004812369
1 0.004812369          NA

Arrived at a reliable solution
```
Note the equal rates in the transition rate matrix: this is what we asked for with the ER model!

Once you've followed the script to check out the format of the ancestral states, plot your results. You should see something like this:

![](https://github.com/joyceem/MPEP_tutorials/blob/90a90e886ffa5f455e03c8b6015a51c2c8228c94/tutorials/images/Pasted%20image%2020250309121027.png)

Now, if we want to get the values for the proportional likelihood (~probability) of each ancestral state, we need to identify the node number, and extract the values for that node.

```
melioideae <- getMRCA(tree, tip=c("MELI_Quivisianthe_papinae", "MELI_Pterorhachis_zenkeri"))

melioideae #This is a node/edge number in ape's language, to get the corresponding row in rayDISC's ans$states, we need to substract the number of tips

nrtips <- length(tree$tip.label)
ans$states[melioideae-nrtips,]
```
What is the most likely seed type in the ancestor Melioideae according to this model? 

Can you do the same for the ancesral node of the Cedreloideae subfamily?

Now, we are going to model evolution again according to a Markov Model, but this time we will use an ARD transition rate matrix.

```
ans <- rayDISC(tree, data, ntraits=1, charnum=1, model="ARD")
```

Note that this time, we have unequal rates in the transition rate matrix (this is what we asked for with the ARD model):
```
> ans
Fit
      -lnL      AIC     AICc ntax
 -10.36366 24.72732 25.22732   27

Rates
            0           1
0          NA 0.002967107
1 0.008803931          NA

Arrived at a reliable solution
```
What does this suggest about how seeds evolve in Meliaceae?

Plot the results of the ARD model.

![](https://github.com/joyceem/MPEP_tutorials/blob/90a90e886ffa5f455e03c8b6015a51c2c8228c94/tutorials/images/Pasted%20image%2020250309122019.png)

How do the results differ from the ER model?

But which one is a better model for our data? We can test this by comparing the AIC:
```
> ansER$AIC
[1] 24.11076
> ansARD$AIC
[1] 24.72732
> deltaAIC <- ansER$AIC - ansARD$AIC
> deltaAIC
[1] -0.6165636
```
We can see that the ER is a slightly better fit than the ARD model.

#### Modelling sexual system evolution (multistate trait)
Follow the R script to repeat what we just did for seeds to reconstruct the ancestral state of sexual systems in Meliaceae.

When you plot your results, you should see something like this for the ER...

![](https://github.com/joyceem/MPEP_tutorials/blob/90a90e886ffa5f455e03c8b6015a51c2c8228c94/tutorials/images/Pasted%20image%2020250309124021.png)

...and this for the ARD:

![](https://github.com/joyceem/MPEP_tutorials/blob/90a90e886ffa5f455e03c8b6015a51c2c8228c94/tutorials/images/Pasted%20image%2020250309124131.png)

Explore the output. Which is the best model? What are the most likely ancestral states for Melioideae, Cedreloideae and Meliaceae?

### Evolution of discrete characters with rate changes ("Hidden Markov Models")
One potential problem with traditional Mk models such as the ones we just ran is that they assume that transition rates are fixed across the whole tree. One way we can allow for rates to change across the tree is by modelling states that allow for different rate categories across the tree. Let's go back to our seed data to try this.

First, let's rerun an ER model with one rate category:
```
ansER_HMM1 <-corHMM(phy = tree, data = seeds_data, rate.cat = 1, model="ER")
ansER_HMM1
```
Inspect the output. You should see something very similar to out previous ER model, as it is the same model:
```
> ansER_HMM1
Fit
      -lnL      AIC     AICc Rate.cat ntax
 -11.05538 24.11076 24.27076        1   27

Legend
  1   2 
"0" "1" 

Rates
            (1,R1)      (2,R1)
(1,R1)          NA 0.004812382
(2,R1) 0.004812382          NA

Arrived at a reliable solution
```

Now, allow for two rate categories (fast=R1 and slow=R2).
```
ansER_HMM2 <-corHMM(phy = tree, data = seeds_data, rate.cat = 2, model="ER")
ansER_HMM2
```
Inspect the output. We can see that even though we have allowed for two rate categories, most transitions are still happening in the same rate class:
```
> ansER_HMM2$solution
             (1,R1)       (2,R1)       (1,R2)       (2,R2)
(1,R1)           NA  0.004812301 2.194261e-09           NA
(2,R1)  0.004812301           NA           NA 2.194261e-09
(1,R2) 78.431733784           NA           NA 4.249587e-06
(2,R2)           NA 78.431733784 4.249587e-06           NA
```

Now repeat for three rate categories.

Compare the AIC of all models. Which one is the most likely?
```
> ansER_HMM1$AIC
[1] 24.11076
> ansER_HMM2$AIC
[1] 30.11076
> ansER_HMM3$AIC
[1] 39.31278
```

We can see that even though we have allowed for the rate categories to differ across the tree, the model with only one category is still the best. This suggests that the rate of seed evolution is constant across the whole tree; but this might also be due to having such a small phylogeny. The results could be different if we had more tips. What results might we expect to see if Cedreloideae had a faster rate of seed evolution than Meiloideae?

## Evolution of discrete characters with state-dependent diversification model (BiSSE)
Now, we are going to use a state-dependent diversification model to model trait evolution. SSE models are powerful approaches for testing the association of a character with diversification rate heterogeneity. You can implement these SSE models in a number of ways (such as RevBayes). We are going to implement a BiSSE model in R using the R package `diversitree`. For more information about the package see the manual: https://cran.r-project.org/web/packages/diversitree/diversitree.pdf.

The Binary State-dependent Speciation and Extinction model is used when we hypothesise that diversification rates could be dependent on a binary trait. For example, we might hypothesise that winged seeds are associated with greater diversification than unwinged seeds. The BiSSE model has six parameters: 
 - q01 and q10, which describe the rate at which the binary character transitions from one state to another, and
 - λ0, μ0, λ1, and μ1 which describe the diversification dynamics while a lineage is in state 0 and 1, respectively.

Given a phylogeny and character values, we can estimate these parameters in R.

![BiSSE](https://github.com/user-attachments/assets/f868b475-5700-4185-b851-40b13a2529ac)

A schematic overview of the BiSSE model from [https://revbayes.github.io/tutorials/sse/bisse-intro.html#bisse_theory](https://revbayes.github.io/tutorials/sse/bisse-intro.html#bisse_theory).

Follow the R code to format the seed trait data and tree. Set up the BiSSE model. You should get an output that looks something like this:

```
BiSSE likelihood function:
  * Parameter vector takes 6 elements:
     - lambda0, lambda1, mu0, mu1, q01, q10
  * Function takes arguments (with defaults)
     - pars: Parameter vector
     - condition.surv [TRUE]: Condition likelihood on survial?
     - root [ROOT.OBS]: Type of root treatment
     - root.p [NULL]: Vector of root state probabilities
     - intermediates [FALSE]: Also return intermediate values?
  * Phylogeny with 27 tips and 26 nodes
     - Taxa: MELI_Aglaia_spectabilis, MELI_Aphanamixis_polystachya, MELI_Cabralea_canjerana, ...
  * References:
     - Maddison et al. (2007) doi:10.1080/10635150701607033
     - FitzJohn et al. (2009) doi:10.1093/sysbio/syp067
R definition:
function (pars, condition.surv = TRUE, root = ROOT.OBS, root.p = NULL, intermediates = FALSE)
```

This gives us information about the setup of our model. You can see that there are the 6 parameters that we will need to optimise to find the best fit for our data.

We could do this manually:

```
> # Compute the likelihood under different parameter values
> ## Parameter order: lambda0, lambda1, mu0, mu1, q01, q10 
> bisse_model(c(1.0, 4.0, 0.5, 0.5, 0.2, 0.8))
[1] -731.7963
> bisse_model(c(4.0, 1.0, 0  , 1  , 0.5, 0.5))
[1] -774.51
```

But this will take forever! So, instead we can use `find.mle()` to find the maximum likelihood estimate (MLE) for our parameters.

```
initial_pars<-starting.point.bisse(tree) # provide starting values
fit_model<-find.mle(bisse_model,initial_pars) # run search algorithm
```
What is the likelihood of the best model? 

```
> fit_model$lnLik
[1] -126.7945
```

Much better than the likelihood of the models when we put in the parameter values by hand!

Now, we can look at what values the 6 parameters had in the optimal model. 

```
> round(fit_model$par,digits=3) ##we round to two digits for easier reading
lambda0 lambda1     mu0     mu1     q01     q10 
  0.024   0.025   0.000   0.000   0.000   0.010
```

What does this tell us about the association between each trait state, and speciation and extinction?

## Evolution of discrete characters with state-dependent diversification models with hidden states (HiSSE)

Now this is ok, but one issue with the BiSSE model is that it has trouble distinguishing whether changes in diversification rates are due to the trait state, or whether changes in diversification rates are independent of that trait and due to something that we haven't observed or captured in our data (i.e. the 'hidden state'). The HiSSE model (hidden-state dependent speciation and extinction) incorporates a 2nd, unobserved trait to account for correlations that are not tied to the character we are testing. Any changes in the hidden trait's state indicate that diversification rate shifts are not correlated with the observed character (rather, they are correlated with something we haven't observed) [[Beaulieu & O’Meara 2016](https://doi.org/10.1093/sysbio/syw022)].

![HiSSE](https://github.com/user-attachments/assets/a9fc99ba-64eb-488a-90ae-efa8fe5d756e)

A schematic overview of the HiSSE model from [https://revbayes.github.io/tutorials/sse/hisse.html](https://revbayes.github.io/tutorials/sse/hisse.html)

The first thing we need to do is set up the transition rate matrix. This is separate from the diversification rate parameters, which reinforces the idea that SSE models are not simply trait evolution models, but rather joint models for the tree *and* the evolution of a trait.

```
> trans.rate.hisse <- TransMatMakerHiSSE(hidden.traits=1)
> print(trans.rate.hisse)
     (0A) (1A) (0B) (1B)
(0A)   NA    2    5   NA
(1A)    1   NA   NA    5
(0B)    5   NA   NA    4
(1B)   NA    5    3   NA

```
The observed characters are denoted with either a 0 or a 1 while the hidden characters are denoted with A and B. The numbers in the matrix (1-5) don’t tell us what the actual rates are, just that elements with the same number will share a rate in our model. Here, we have set it so that  transitions from 1→0 all have the same rate,  and transitions from 0→1 have the same rate. We can also see a few NAs in our matrix: we see these whenever we see both states transition at once ( 1A↔0B or 0A↔1B). Generally, we assume that only one state changes at a given time so we want to exclude the possibility that both transitions occur at once.

Next we need to specify which diversification rates are shared between states. Unlike `diversitree`, `hisse` has uses a different parameterization for diversification. Instead of speciation and extinction, `hisse` uses net turnover and the extinction fraction. These are just transformations of speciation and extinction, and in effect can be interpreted as speciation and extinction rates. Specifically they are:
- net turnover=λ+μ
- extinction fraction=μ/λ

The states are ordered 0A, 1A, 0B, 1B, so we will specify that the diversification rates are different between 0 and 1, and can differ between A and B - i.e. that all states 0A and 1A are different.

```
turnover <- c(1,2,3,4)
extinction.fraction <- c(1,2,3,4)
```

Now we have set up how all of the parameters for the model can behave, we can run it:
```
HiSSE_model <- hisse(phy=tree, data=seeds_data, turnover=turnover, 
               eps=extinction.fraction, hidden.states=TRUE, 
               trans.rate=trans.rate.hisse)
```

Check out the results. What do they show?

```
> HiSSE_model
Fit 
            lnL             AIC            AICc          n.taxa n.hidden.states 
      -123.3304        272.6607        300.6607         27.0000          2.0000 

Model parameters: 

  turnover0A   turnover1A        eps0A        eps1A        q0A1A        q1A0A        q0A0B        q1A1B   turnover0B   turnover1B 
2.061154e-09 2.923958e-02 3.000000e+00 2.869615e-09 2.061154e-09 2.061154e-09 6.064881e-02 6.064881e-02 4.231171e-02 1.723908e-02 
       eps0B        eps1B        q0B1B        q1B0B        q0B0A        q1B1A 
2.061154e-09 2.061154e-09 2.061154e-09 2.416917e-02 6.064881e-02 6.064881e-02
```

Now we can actually reconstruct the likelihood of all of the ancestral states and plot them to see what it looks like on our tree. 

![](https://github.com/joyceem/MPEP_tutorials/blob/90a90e886ffa5f455e03c8b6015a51c2c8228c94/tutorials/images/Pasted%20image%2020250309135744.png)

Now, think about these results carefully. Do they actually make sense? Just because a model is 'fancier' because it has more parameters, and even if it has a high likelihood statistically, you should always carefully consider whether the model and its results actually _make sense_. I would suggest that in this case we have an example of overparametisation for the data we have. We have a relatively small dataset, with a complex model. There is probably not enough tips in this tree to come up with a meaningful answer. Also note that we have only sampled at the genus level, not the species level; if we included species we would be much better placed to estimate diversification because we have species as our units, not genera (which are much more arbitrary, and might not accurately represent extant diversity).
