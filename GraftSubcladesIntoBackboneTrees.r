# =============================================================================
# GraftSubcladesIntoBackboneTrees.r
#
# Author:  Ricardo Betancur (adapted for general use)
# Date:    March 2026
# Purpose: Graft one or more densely-sampled subclade phylogenies into one or
#          more backbone (spine) trees, producing every valid combination of
#          grafted trees. Useful for integrating clade-specific phylogenies
#          into a set of backbone trees.
#
# HOW IT WORKS (the short version)
# ---------------------------------
# 1. You provide backbone tree(s) and subclade tree(s) in NEXUS or Newick format.
# 2. You specify two "key taxa" — tips that appear in BOTH the backbone and the
#    subclade. These anchor the graft:
#      - Key.taxon1 marks WHERE on the backbone the subclade will attach.
#      - Key.taxon2 (together with Key.taxon1) defines the MRCA in the subclade
#        so we can extract the clade of interest.
# 3. The script:
#      a. Finds the MRCA of Key.taxon1 + Key.taxon2 in the subclade tree and
#         extracts that clade.
#      b. Drops Key.taxon2 from the backbone (so only Key.taxon1 remains as the
#         attachment point).
#      c. Uses bind.tree() to attach the extracted subclade at Key.taxon1's
#         position, placing it at the crown age of the subclade.
#      d. Removes Key.taxon1 (the placeholder tip) from the grafted result.
# 4. This is repeated for every backbone x subclade combination.
# 5. All grafted trees are saved to a single output file.
#
# DEPENDENCIES
# ------------
# install.packages("phytools")  # if not already installed
# install.packages("ape")       # if not already installed
#
# =============================================================================


# ---- Load required packages -------------------------------------------------

library(ape)       # core phylogenetics: read/write trees, drop.tip, bind.tree
library(phytools)  # fastMRCA, extract.clade, node.depth.edgelength, etc.


# =============================================================================
#                        USER CONFIGURATION SECTION
#
#  Edit the values below to match YOUR files and taxa, then source this script.
# =============================================================================

# -- 1. Input files -----------------------------------------------------------
# Paths to your tree files. Each file can contain one or many trees.
# Accepted formats: NEXUS (.nex, .nex.tre) or Newick (.tre, .nwk).

# NOTE: These use if(!exists(...)) so you can set them BEFORE source()-ing this
# script (e.g., from a wrapper like run_example.r) and they won't be overwritten.

if (!exists("BACKBONE_FILE")) BACKBONE_FILE <- "backbone_trees.nex"      # <-- your backbone tree file
if (!exists("SUBCLADE_FILE")) SUBCLADE_FILE <- "subclade_trees.nex"       # <-- your subclade tree file

# -- 2. Output file -----------------------------------------------------------

if (!exists("OUTPUT_FILE")) OUTPUT_FILE <- "Grafted_Trees.nex"          # <-- name for output file (NEXUS)

# -- 3. Key taxa --------------------------------------------------------------
# Two taxa that are present in BOTH the backbone and the subclade trees.
# They define the grafting point:
#   Key.taxon1 = the tip in the backbone where the subclade will be attached.
#   Key.taxon2 = used together with Key.taxon1 to find the MRCA in the subclade.
#                Also dropped from the backbone before grafting.
#
# Example: if grafting a family-level phylogeny into an order-level backbone,
# pick two species from that family that appear in both trees.

if (!exists("Key.taxon1")) Key.taxon1 <- "Genus_species1"             # <-- REPLACE with your taxon
if (!exists("Key.taxon2")) Key.taxon2 <- "Genus_species2"             # <-- REPLACE with your taxon

# -- 4. Optional: rescale subclade branch lengths? ----------------------------
# If your subclade branch lengths are on a different scale than the backbone
# (e.g., substitutions x100 vs. time-calibrated), set RESCALE to TRUE and
# choose a divisor. Otherwise leave as FALSE.

if (!exists("RESCALE_SUBCLADE")) RESCALE_SUBCLADE <- FALSE
if (!exists("RESCALE_DIVISOR"))  RESCALE_DIVISOR  <- 100                    # branch lengths will be divided by this


# =============================================================================
#                           FUNCTION DEFINITIONS
# =============================================================================

# -----------------------------------------------------------------------------
# read_trees: Load trees from a file, auto-detecting NEXUS vs. Newick format.
#
# Arguments:
#   filepath  Path to tree file (.nex, .nex.tre, .tre, .nwk)
#
# Returns:
#   A multiPhylo object (list of trees), even if the file has just one tree.
# -----------------------------------------------------------------------------
read_trees <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(paste("File not found:", filepath))
  }

  # Peek at first line to decide format
  first_line <- toupper(trimws(readLines(filepath, n = 1)))

  if (grepl("^#NEXUS", first_line)) {
    trees <- read.nexus(filepath)
  } else {
    trees <- read.tree(filepath)
  }

  # Ensure multiPhylo even if single tree
  if (inherits(trees, "phylo")) {
    trees <- c(trees)                    # converts phylo -> multiPhylo
    class(trees) <- "multiPhylo"
    names(trees) <- tools::file_path_sans_ext(basename(filepath))
  }

  cat(sprintf("  Loaded %d tree(s) from '%s'\n", length(trees), filepath))
  return(trees)
}


# -----------------------------------------------------------------------------
# rescale_branch_lengths: Divide all branch lengths in a set of trees by a
#   constant. Useful when subclade trees are in substitutions (x100) and the
#   backbone is time-calibrated.
#
# Arguments:
#   trees    multiPhylo object
#   divisor  number to divide branch lengths by (default = 100)
#
# Returns:
#   A rescaled multiPhylo object.
# -----------------------------------------------------------------------------
rescale_branch_lengths <- function(trees, divisor = 100) {
  rescaled <- lapply(trees, function(tree) {
    tree$edge.length <- tree$edge.length / divisor
    return(tree)
  })
  class(rescaled) <- "multiPhylo"
  names(rescaled) <- names(trees)
  return(rescaled)
}


# -----------------------------------------------------------------------------
# graft_all_combinations: The core grafting engine.
#
# For every subclade tree and every backbone tree, this function:
#   1. Extracts the subclade of interest using the MRCA of Key.taxon1 + Key.taxon2.
#   2. Finds the MRCA of those two key taxa in the backbone and drops ALL
#      descendant tips EXCEPT Key.taxon1. This collapses the clade down to a
#      single placeholder tip whose pendant edge stretches back to the MRCA node,
#      giving bind.tree() enough branch length to place the subclade root at
#      the correct depth.
#   3. Attaches the extracted subclade at Key.taxon1 using bind.tree().
#   4. Removes the placeholder tip (Key.taxon1) from the final tree.
#
# Arguments:
#   Backbone.trees   multiPhylo: one or more backbone trees
#   Subclade.trees   multiPhylo: one or more subclade trees
#   Key.taxon1       character: tip present in both (attachment point)
#   Key.taxon2       character: tip present in both (defines MRCA in subclade;
#                    removed from backbone before grafting)
#
# Returns:
#   A named list of grafted phylo objects. Names follow the pattern:
#   "BackboneName_SubcladeName" so you can trace which inputs produced each tree.
# -----------------------------------------------------------------------------
graft_all_combinations <- function(Backbone.trees, Subclade.trees,
                                   Key.taxon1, Key.taxon2) {

  grafted_trees <- list()

  # -- Inner helper: performs the actual bind.tree + cleanup ------------------
  graft_one <- function(Backbone, Subtree_clade, Keynode_name) {
    # --- Calculate the correct graft position ---
    #
    # The position parameter in bind.tree() is measured FROM the key taxon's
    # tip going rootward along its pendant edge. It determines where the
    # subclade root will be placed.
    #
    # For tip-dated trees, Key.taxon1 may be a fossil (not at the present).
    # The correct position must ensure that EXTANT tips in the grafted
    # subclade align at the present (t = 0) with extant tips in the backbone.
    #
    # Formula:
    #   position = subclade_root_age - key_taxon_age_in_backbone
    #
    # Where:
    #   subclade_root_age     = max root-to-tip distance in the subclade
    #                           (= age of subclade root from the present)
    #   key_taxon_age_in_backbone = tree_height - depth_of_key_taxon
    #                           (= how far the key taxon is from the present)
    #
    # For ultrametric trees this simplifies to the crown age, since all tips
    # (including Key.taxon1) are at the present (age = 0).

    # Subclade root age (from the present = from extant tips)
    sub_depths        <- node.depth.edgelength(Subtree_clade)
    subclade_root_age <- max(sub_depths)  # root is at depth 0, extant tips at max

    # Key taxon age in the backbone (how old is this tip?)
    bb_depths         <- node.depth.edgelength(Backbone)
    bb_height         <- max(bb_depths)
    keynode_idx       <- which(Backbone$tip.label == Keynode_name)
    key_age_bb        <- bb_height - bb_depths[keynode_idx]

    # Position from the key taxon tip going rootward
    graft_pos <- subclade_root_age - key_age_bb

    # Safety check: position must be positive and fit within the branch
    edge_row   <- which(Backbone$edge[, 2] == keynode_idx)
    branch_len <- Backbone$edge.length[edge_row]

    if (graft_pos <= 0) {
      warning(sprintf("Subclade root age (%.2f) <= key taxon age in backbone (%.2f). Using small positive position.",
                      subclade_root_age, key_age_bb))
      graft_pos <- branch_len * 0.01
    } else if (graft_pos > branch_len) {
      cat(sprintf("    Note: needed position (%.2f) > branch length (%.2f); capping.\n",
                  graft_pos, branch_len))
      graft_pos <- branch_len
    }

    # Attach the subclade at that tip position
    grafted <- bind.tree(Backbone, Subtree_clade,
                         where = keynode_idx,
                         position = graft_pos,
                         interactive = FALSE)

    # Remove the placeholder tip — its role is done
    grafted <- drop.tip(grafted, Keynode_name)

    return(grafted)
  }

  # -- Loop over every subclade x backbone combination -----------------------
  for (s in seq_along(Subclade.trees)) {
    sub_tree   <- Subclade.trees[[s]]
    sub_label  <- names(Subclade.trees)[s]
    if (is.null(sub_label)) sub_label <- paste0("Subclade_", s)

    # Validate that key taxa exist in this subclade tree
    if (!all(c(Key.taxon1, Key.taxon2) %in% sub_tree$tip.label)) {
      warning(sprintf("Skipping subclade '%s': key taxa not both present.", sub_label))
      next
    }

    # Find the MRCA of the two key taxa in the subclade tree, then extract
    # that clade. This is the piece we will graft into the backbone.
    mrca_node     <- fastMRCA(sub_tree, Key.taxon1, Key.taxon2)
    subclade_part <- extract.clade(sub_tree, mrca_node)

    for (b in seq_along(Backbone.trees)) {
      bb_tree  <- Backbone.trees[[b]]
      bb_label <- names(Backbone.trees)[b]
      if (is.null(bb_label)) bb_label <- paste0("Backbone_", b)

      # Validate that key taxa exist in this backbone tree
      if (!all(c(Key.taxon1, Key.taxon2) %in% bb_tree$tip.label)) {
        warning(sprintf("Skipping backbone '%s': key taxa not both present.", bb_label))
        next
      }

      # --- Critical step: prune ALL subclade members from backbone except
      #     Key.taxon1. This ensures Key.taxon1's pendant edge extends all
      #     the way back to the MRCA node, so bind.tree() has enough branch
      #     length for position = crown_age. ---
      mrca_bb       <- fastMRCA(bb_tree, Key.taxon1, Key.taxon2)
      clade_tips_bb <- extract.clade(bb_tree, mrca_bb)$tip.label
      tips_to_drop  <- setdiff(clade_tips_bb, Key.taxon1)
      pruned_bb     <- drop.tip(bb_tree, tips_to_drop)

      # Graft!
      result <- tryCatch({
        graft_one(pruned_bb, subclade_part, Key.taxon1)
      }, error = function(e) {
        warning(sprintf("Grafting failed for '%s' x '%s': %s",
                        bb_label, sub_label, e$message))
        return(NULL)
      })

      if (!is.null(result)) {
        combo_name <- paste0(bb_label, "_", sub_label)
        grafted_trees[[combo_name]] <- result
      }
    }
  }

  cat(sprintf("  Produced %d grafted tree(s).\n", length(grafted_trees)))
  return(grafted_trees)
}


# =============================================================================
#                              MAIN WORKFLOW
# =============================================================================

cat("\n========== Grafting Subclades into Backbone Trees ==========\n\n")

# Step 1: Load trees
cat("Step 1: Loading trees...\n")
Backbone.trees <- read_trees(BACKBONE_FILE)
Subclade.trees <- read_trees(SUBCLADE_FILE)

# Step 2: (Optional) Rescale subclade branch lengths
if (RESCALE_SUBCLADE) {
  cat(sprintf("Step 2: Rescaling subclade branch lengths (dividing by %d)...\n",
              RESCALE_DIVISOR))
  Subclade.trees <- rescale_branch_lengths(Subclade.trees, RESCALE_DIVISOR)
} else {
  cat("Step 2: Skipping rescaling (RESCALE_SUBCLADE = FALSE).\n")
}

# Step 3: Validate key taxa are present
cat("Step 3: Checking key taxa...\n")
# Quick check on the first tree of each set
bb_tips  <- Backbone.trees[[1]]$tip.label
sub_tips <- Subclade.trees[[1]]$tip.label

if (!(Key.taxon1 %in% bb_tips))
  stop(paste("Key.taxon1 not found in backbone:", Key.taxon1))
if (!(Key.taxon2 %in% bb_tips))
  stop(paste("Key.taxon2 not found in backbone:", Key.taxon2))
if (!(Key.taxon1 %in% sub_tips))
  stop(paste("Key.taxon1 not found in subclade:", Key.taxon1))
if (!(Key.taxon2 %in% sub_tips))
  stop(paste("Key.taxon2 not found in subclade:", Key.taxon2))
cat("  Key taxa verified in both tree sets.\n")

# Step 4: Graft all combinations
cat("Step 4: Grafting all backbone x subclade combinations...\n")
all_grafted <- graft_all_combinations(Backbone.trees, Subclade.trees,
                                      Key.taxon1, Key.taxon2)

# Step 5: Save results
if (length(all_grafted) > 0) {
  cat(sprintf("Step 5: Writing %d grafted trees to '%s'...\n",
              length(all_grafted), OUTPUT_FILE))
  write.nexus(all_grafted, file = OUTPUT_FILE)
  cat(sprintf("\nDone! %d grafted trees saved to: %s\n", length(all_grafted), OUTPUT_FILE))
} else {
  cat("\nNo trees were successfully grafted. Check warnings above.\n")
}


# =============================================================================
#                     QUICK VISUAL CHECK (optional)
# =============================================================================
# Uncomment the lines below to plot the first grafted tree as a sanity check.
#
# if (length(all_grafted) > 0) {
#   plot(all_grafted[[1]], cex = 0.4, no.margin = TRUE)
#   title(main = names(all_grafted)[1], cex.main = 0.8)
# }


# =============================================================================
#                          TIPS & TROUBLESHOOTING
# =============================================================================
#
# Q: How do I pick Key.taxon1 and Key.taxon2?
# A: They must be tips present in BOTH the backbone and the subclade tree.
#    Ideally, they span the root of the subclade you want to graft — i.e., they
#    are the two most distantly-related taxa in the subclade that also appear in
#    the backbone. Their MRCA in the subclade defines which part gets extracted.
#
# Q: My subclade tree contains outgroups I don't want grafted.
# A: That's fine! The script uses extract.clade() on the MRCA of the two key
#    taxa, so only the ingroup (defined by those taxa) gets grafted. Anything
#    outside the MRCA is ignored.
#
# Q: What if I have multiple subclades to graft sequentially?
# A: Run this script once for each subclade. After the first graft, use the
#    OUTPUT_FILE as the new BACKBONE_FILE for the next round. Example:
#
#      Round 1: Backbone.nex + Spikefishes.nex  -> Grafted_round1.nex
#      Round 2: Grafted_round1.nex + OtherClade.nex -> Grafted_round2.nex
#
# Q: What if branch lengths don't match (e.g., substitutions vs. time)?
# A: Set RESCALE_SUBCLADE <- TRUE and adjust RESCALE_DIVISOR so the subclade
#    branch lengths are on the same scale as the backbone.
#
# Q: I get "Key taxa not both present" warnings.
# A: Double-check the exact spelling of your key taxa (tip labels are
#    case-sensitive and often use underscores instead of spaces).
#    Run: Backbone.trees[[1]]$tip.label  to see all tip names.
#
# =============================================================================
