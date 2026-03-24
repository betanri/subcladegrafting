# Graft Subclades into Backbone Phylogenies

R script for grafting densely-sampled subclade phylogenies into larger backbone trees, producing all valid combinations of grafted trees.

## What it does

Given a set of backbone trees and a set of subclade trees that share two "key taxa," the script:

1. Finds the MRCA of the two key taxa in each subclade tree and extracts that clade
2. Removes one key taxon from the backbone, leaving the other as the attachment point
3. Attaches the extracted subclade at the remaining key taxon using `bind.tree()`, placing the junction at the subclade's crown age
4. Removes the placeholder tip from the final tree
5. Repeats for every backbone x subclade combination

This is useful for integrating clade-specific phylogenies (e.g., from targeted-capture or UCE datasets) into a broader reference tree (e.g., a time-calibrated supermatrix).

## Requirements

```r
install.packages("ape")
install.packages("phytools")
```

## Quick start

1. Clone this repo
2. Edit four lines at the top of `GraftSubcladesIntoBackboneTrees.r`:

```r
BACKBONE_FILE <- "my_backbone.nex"       # your backbone tree file (NEXUS or Newick)
SUBCLADE_FILE <- "my_subclade.nex"        # your subclade tree file
Key.taxon1    <- "Genus_species1"         # tip in BOTH trees (attachment point)
Key.taxon2    <- "Genus_species2"         # tip in BOTH trees (defines MRCA; dropped from backbone)
```

3. Source the script:

```r
source("GraftSubcladesIntoBackboneTrees.r")
```

Output is written as a NEXUS file containing all grafted trees.

## Example

The `example/` directory includes backbone and subclade trees (Liparidae snailfishes within a percomorph backbone) with a pre-configured runner:

```r
setwd("example")
source("run_example.r")
```

This grafts 10 subclade trees into 10 backbone trees, producing 100 grafted trees.

## Choosing key taxa

The two key taxa must be tip labels present in **both** the backbone and subclade trees. Pick two species that span the root of the subclade you want to graft (i.e., distantly-related within the clade). Their MRCA in the subclade tree defines which portion gets extracted.

To see shared tips between your trees:

```r
library(ape)
bb  <- read.nexus("backbone.nex")
sub <- read.nexus("subclade.nex")
intersect(bb[[1]]$tip.label, sub[[1]]$tip.label)
```

## Grafting multiple subclades sequentially

Run the script once per subclade, using the previous output as the new backbone:

```
Round 1:  Backbone.nex   + SubcladeA.nex  ->  Grafted_round1.nex
Round 2:  Grafted_round1.nex + SubcladeB.nex  ->  Grafted_round2.nex
```

## Branch-length rescaling

If the subclade branch lengths are on a different scale than the backbone (e.g., substitutions vs. time), set:

```r
RESCALE_SUBCLADE <- TRUE
RESCALE_DIVISOR  <- 100   # branch lengths divided by this value
```

## Citation

If you use this script, please cite:

Betancur-R., R. et al. (in prep). Phylogenetic comparative methods for Australian and New Guinean freshwater fishes.
