# run_example.r
#
# Demonstrates grafting using the included example data (Liparidae snailfishes
# grafted into a percomorph backbone). Run from the example/ directory:
#
#   setwd("example")
#   source("run_example.r")
#
# The backbone contains 10 trees spanning several percomorph orders, each with
# a sparse Liparidae sampling. The subclade file contains 10 densely-sampled
# Liparidae trees. The script grafts every backbone x subclade combination
# (10 x 10 = 100 grafted trees).
# -----------------------------------------------------------------------------

# ---- Configuration (already filled in for the example data) -----------------

BACKBONE_FILE <- "backbone_trees.nex"
SUBCLADE_FILE <- "subclade_trees.nex"
OUTPUT_FILE   <- "Example_Grafted_Trees.nex"

# Two Liparidae species present in both backbone and subclade trees.
# They span the root of Liparidae, so extract.clade() captures the full family.
Key.taxon1 <- "Eupercaria_Perciformes_Liparidae_Careproctus_gilberti_target"
Key.taxon2 <- "Eupercaria_Perciformes_Liparidae_Liparis_gibbus_target"

RESCALE_SUBCLADE <- FALSE

# ---- Source the main grafting script (one directory up) ---------------------

source("../GraftSubcladesIntoBackboneTrees.r")

# ---- Quick sanity check: plot the first grafted tree -----------------------

if (length(all_grafted) > 0) {
  cat("\nPlotting first grafted tree as a sanity check...\n")
  plot(all_grafted[[1]], cex = 0.3, no.margin = TRUE)
  title(main = names(all_grafted)[1], cex.main = 0.7)
}
