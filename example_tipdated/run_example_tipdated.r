# run_example_tipdated.r
#
# Demonstrates grafting with TIP-DATED (non-ultrametric) trees.
# These are time-calibrated trees where extinct/fossil taxa have tips
# at different time points, so they are NOT ultrametric by design.
#
# Run from the example_tipdated/ directory:
#
#   setwd("example_tipdated")
#   source("run_example_tipdated.r")
#
# The backbone is a single Tetraodontiformes MCC tree (246 tips).
# The subclade is a single spikefish (Triacanthodidae) tree (26 tips).
# Both include fossil taxa with non-zero tip ages.
# -----------------------------------------------------------------------------

# ---- Configuration ---------------------------------------------------------

BACKBONE_FILE <- "backbone_tipdated.nwk"
SUBCLADE_FILE <- "subclade_tipdated.nwk"
OUTPUT_FILE   <- "Grafted_TipDated.nex"

# Two species present in both trees that span the root of the subclade.
# Key.taxon1 is a fossil (~24 Ma) — the script correctly handles fossil
# key taxa by adjusting the graft position so extant tips align at the present.
Key.taxon1 <- "Prohollardia_avita"
Key.taxon2 <- "Atrophacanthus_japonicus"

RESCALE_SUBCLADE <- FALSE

# ---- Source the main grafting script (one directory up) ---------------------

source("../GraftSubcladesIntoBackboneTrees.r")

# ---- Verification ----------------------------------------------------------

cat("\n========== TIP-DATED TREE CHECK ==========\n")
tree <- all_grafted[[1]]
cat(sprintf("Grafted tree tips: %d\n", Ntip(tree)))
cat(sprintf("Is ultrametric: %s  (expected FALSE for tip-dated trees)\n",
            is.ultrametric(tree)))
cat(sprintf("Has branch lengths: %s\n", !is.null(tree$edge.length)))

# Verify extant tips in subclade align with extant tips in backbone
g_depths <- node.depth.edgelength(tree)
g_max    <- max(g_depths)
sub_extant_depth <- g_depths[which(tree$tip.label == "Triacanthodes_anomalus")]
bb_extant_depth  <- g_depths[which(tree$tip.label == "Balistes_capriscus")]
cat(sprintf("Extant tip alignment (subclade vs backbone): %.6f vs %.6f  (diff=%.2e)\n",
            sub_extant_depth, bb_extant_depth, abs(sub_extant_depth - bb_extant_depth)))

# Root-to-tip distance range (will vary for tip-dated trees)
depths <- node.depth.edgelength(tree)
tip_depths <- depths[1:Ntip(tree)]
cat(sprintf("Root-to-tip range: %.2f -- %.2f\n", min(tip_depths), max(tip_depths)))

# ---- Quick sanity check: plot the first grafted tree -----------------------

if (length(all_grafted) > 0) {
  cat("\nPlotting grafted tree as a sanity check...\n")
  plot(all_grafted[[1]], cex = 0.3, no.margin = TRUE)
  title(main = "Tip-dated graft: Triacanthodidae into Tetraodontiformes",
        cex.main = 0.7)
}
