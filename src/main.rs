use rust_python_tree_distances::io::{read_beast_trees};
use std::path::PathBuf;

/// Load a NEXUS file, extract all lines beginning with `tree` (case-insensitive),
/// parse the Newick portion after the '=' and build a vector of phylogenetic trees.
fn main() {
    // Path to NEXUS file with multiple tree definitions.
    let path = PathBuf::from("tests/data/hiv1.trees");
    let burnin_trees = 0;
    let burnin_states = 0;
    let use_real_taxa = true; // replace numeric taxon IDs with labels when available

    let trees = read_beast_trees(&path, burnin_trees, burnin_states, use_real_taxa);

    println!("Loaded {} trees from {:?}", trees.len(), path);
    // // Compute pairwise Robinson–Foulds distances.
    // let distances = compute_pairwise_robinson_foulds(&trees);
    // println!("Pairwise Robinson–Foulds distances (i, j, d):");
    // for (i, j, d) in distances.iter().take(25) {
    //     println!("{i}\t{j}\t{d}");
    // }
    // if distances.len() > 25 {
    //     println!("... ({} total pairs)", distances.len());
    // }
}
