use std::time::Instant;
use rust_python_tree_distances::distances::{
    compute_pairwise_rf_custom_seq,
    compute_pairwise_rf_custom_parallel,
};
use phylotree::tree::Tree as PhyloTree;

fn make_trees() -> Vec<PhyloTree> {
    // A small suite of moderately complex rooted trees (8 taxa)
    let newicks = [
        "(A:0.1,(B:0.2,(C:0.3,D:0.4):0.5):0.6,(E:0.7,(F:0.8,G:0.9):0.95):1.0,H:1.1)R;",
        "((A:0.1,B:0.2):0.25,(C:0.3,(D:0.4,E:0.5):0.55):0.6,(F:0.8,(G:0.9,H:1.0):0.85):0.9)R;",
        "(A:0.1,B:0.2,(C:0.3,(D:0.4,E:0.5):0.55):0.6,(F:0.8,(G:0.9,H:1.0):0.85):0.9)R;",
        "((A:0.1,(B:0.2,C:0.3):0.35):0.4,(D:0.4,E:0.5,(F:0.8,(G:0.9,H:1.0):0.85):0.9):0.95)R;",
        "((A:0.15,B:0.22):0.3,(C:0.28,(D:0.37,E:0.52):0.6):0.62,(F:0.78,(G:0.88,H:1.02):0.82):0.92)R;",
        "(A:0.12,(B:0.21,(C:0.33,D:0.47):0.51):0.58,(E:0.73,(F:0.84,G:0.91):0.96):1.02,H:1.08)R;",
    ];
    let mut trees: Vec<PhyloTree> = newicks.iter().map(|s| PhyloTree::from_newick(s).unwrap()).collect();
    // Replicate to increase problem size slightly
    let extra: Vec<PhyloTree> = newicks.iter().map(|s| PhyloTree::from_newick(s).unwrap()).collect();
    trees.extend(extra);
    trees
}

fn main() {
    let trees = make_trees();
    let reps = 1usize;

    // Warmup to avoid first-run overheads skewing results
    let _ = compute_pairwise_rf_custom_seq(&trees).unwrap();
    let _ = compute_pairwise_rf_custom_parallel(&trees).unwrap();

    // Sequential
    let t0 = Instant::now();
    for _ in 0..reps {
        let _ = compute_pairwise_rf_custom_seq(&trees).unwrap();
    }
    let seq_dur = t0.elapsed();

    // Parallel
    let t1 = Instant::now();
    for _ in 0..reps {
        let _ = compute_pairwise_rf_custom_parallel(&trees).unwrap();
    }
    let par_dur = t1.elapsed();

    println!("Trees: {}  Pairs: {}  Reps: {}", trees.len(), trees.len() * (trees.len()-1) / 2, reps);
    println!("Sequential total: {:?}", seq_dur);
    println!("Parallel   total: {:?}", par_dur);
    let speedup = (seq_dur.as_secs_f64() / par_dur.as_secs_f64()).max(0.0);
    println!("Speedup: {:.2}x", speedup);
}
