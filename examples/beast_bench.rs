use std::time::Instant;

use phylotree::tree::Tree as PhyloTree;
use rust_python_tree_distances::distances::{
    build_bitset_snapshots,
    compute_pairwise_kf_bitset_parallel,
    compute_pairwise_kf_bitset_seq,
    compute_pairwise_rf_bitset_parallel,
    compute_pairwise_rf_bitset_seq,
    compute_pairwise_weighted_rf_bitset_parallel,
    compute_pairwise_weighted_rf_bitset_seq,
};
use rust_python_tree_distances::io::read_beast_trees;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Metric { Rf, Weighted, Kf }

// ==== Hardcoded configuration (edit these to your environment) ====
const PATH: &str = "tests/data/LASV_GENBANK_20250725_curated.L.noNIG.beast2.trees";
const BURNIN_TREES: usize = 0;
const BURNIN_STATES: usize = 3300000;
const USE_REAL_TAXA: bool = false;
const LIMIT: usize = 0; // 0 = no limit; e.g., 1000 to cap
const REPS: usize = 1;   // times to repeat measurements
const METRIC: Metric = Metric::Rf; // Metric::Rf | Metric::Weighted | Metric::Kf

fn main() {
    println!(
        "Reading BEAST trees: path={} burnin_trees={} burnin_states={} use_real_taxa={} limit={}",
        PATH, BURNIN_TREES, BURNIN_STATES, USE_REAL_TAXA, LIMIT
    );
    let mut trees: Vec<PhyloTree> = read_beast_trees(PATH, BURNIN_TREES, BURNIN_STATES, USE_REAL_TAXA);
    if trees.is_empty() {
        eprintln!("No trees parsed; check the PATH constant and file format.");
        std::process::exit(2);
    }
    if LIMIT > 0 && trees.len() > LIMIT { trees.truncate(LIMIT); }

    let pairs = trees.len() * (trees.len() - 1) / 2;
    println!("Loaded trees: {}  pairs: {}  metric: {:?}  reps: {}", trees.len(), pairs, METRIC, REPS);

    // Build snapshots once (true Rustacean: reuse and test pairwise-only work)
    let snaps = build_bitset_snapshots(&trees).expect("snapshot build");
    match METRIC {
        Metric::Rf => {
            let s1 = compute_pairwise_rf_bitset_seq(&snaps);
            let p1 = compute_pairwise_rf_bitset_parallel(&snaps);
            assert_eq!(s1, p1, "parallel vs sequential RF mismatch");
            let t0 = Instant::now();
            for _ in 0..REPS { let _ = compute_pairwise_rf_bitset_seq(&snaps); }
            let seq_dur = t0.elapsed();
            let t1 = Instant::now();
            for _ in 0..REPS { let _ = compute_pairwise_rf_bitset_parallel(&snaps); }
            let par_dur = t1.elapsed();
            println!("Sequential total: {:?}", seq_dur);
            println!("Parallel   total: {:?}", par_dur);
            println!("Speedup: {:.2}x", seq_dur.as_secs_f64() / par_dur.as_secs_f64());
        }
        Metric::Weighted => {
            let s1 = compute_pairwise_weighted_rf_bitset_seq(&snaps);
            let p1 = compute_pairwise_weighted_rf_bitset_parallel(&snaps);
            assert_eq!(s1.len(), p1.len());
            for (a,b) in s1.iter().zip(p1.iter()) { assert_eq!(a.0, b.0); assert_eq!(a.1, b.1); assert!((a.2 - b.2).abs() < 1e-12); }
            let t0 = Instant::now();
            for _ in 0..REPS { let _ = compute_pairwise_weighted_rf_bitset_seq(&snaps); }
            let seq_dur = t0.elapsed();
            let t1 = Instant::now();
            for _ in 0..REPS { let _ = compute_pairwise_weighted_rf_bitset_parallel(&snaps); }
            let par_dur = t1.elapsed();
            println!("Sequential total: {:?}", seq_dur);
            println!("Parallel   total: {:?}", par_dur);
            println!("Speedup: {:.2}x", seq_dur.as_secs_f64() / par_dur.as_secs_f64());
        }
        Metric::Kf => {
            let s1 = compute_pairwise_kf_bitset_seq(&snaps);
            let p1 = compute_pairwise_kf_bitset_parallel(&snaps);
            assert_eq!(s1.len(), p1.len());
            for (a,b) in s1.iter().zip(p1.iter()) { assert_eq!(a.0, b.0); assert_eq!(a.1, b.1); assert!((a.2 - b.2).abs() < 1e-12); }
            let t0 = Instant::now();
            for _ in 0..REPS { let _ = compute_pairwise_kf_bitset_seq(&snaps); }
            let seq_dur = t0.elapsed();
            let t1 = Instant::now();
            for _ in 0..REPS { let _ = compute_pairwise_kf_bitset_parallel(&snaps); }
            let par_dur = t1.elapsed();
            println!("Sequential total: {:?}", seq_dur);
            println!("Parallel   total: {:?}", par_dur);
            println!("Speedup: {:.2}x", seq_dur.as_secs_f64() / par_dur.as_secs_f64());
        }
    }
}
