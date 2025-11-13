//! Parallel re-implementations of RF, weighted RF, and Kuhner–Felsenstein using
//! only public `phylotree` APIs. Safe to parallelize (no shared Tree state).

use phylotree::tree::Tree as PhyloTree;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

// =============================
// Custom parallel re-implementation (no reliance on private bipartition cache)
// =============================

#[derive(Debug, Clone)]
struct PartitionsSnapshot {
    parts: HashSet<Vec<usize>>,         // canonical partitions as sorted leaf indices
    lengths: HashMap<Vec<usize>, f64>,  // partition -> branch length
    root_children: HashSet<Vec<usize>>, // for rooted RF adjustment
    rooted: bool,
}

fn snapshot_tree_partitions(tree: &PhyloTree) -> Result<PartitionsSnapshot, phylotree::tree::TreeError> {
    use phylotree::tree::TreeError;
    let rooted = tree.is_rooted()?;
    let leaves = tree.get_leaves();
    if leaves.is_empty() { return Err(TreeError::IsEmpty); }

    // Map leaf node id -> compact index [0..n)
    let mut leaf_pos: HashMap<usize, usize> = HashMap::with_capacity(leaves.len());
    for (i, id) in leaves.iter().enumerate() { leaf_pos.insert(*id, i); }

    // DFS to collect leaf sets per node
    fn collect(node_id: usize, tree: &PhyloTree, leaf_pos: &HashMap<usize, usize>, cache: &mut HashMap<usize, Vec<usize>>) -> Vec<usize> {
        if let Some(v) = cache.get(&node_id) { return v.clone(); }
        let node = tree.get(&node_id).expect("valid node");
        if node.children.is_empty() {
            let mut v = Vec::with_capacity(1);
            v.push(*leaf_pos.get(&node_id).expect("leaf mapped"));
            cache.insert(node_id, v.clone());
            return v;
        }
        let mut agg: Vec<usize> = Vec::new();
        for &ch in &node.children { let mut sub = collect(ch, tree, leaf_pos, cache); agg.append(&mut sub); }
        agg.sort_unstable();
        cache.insert(node_id, agg.clone());
        agg
    }

    let mut cache: HashMap<usize, Vec<usize>> = HashMap::new();
    let root_id = tree.get_root()?;
    collect(root_id, tree, &leaf_pos, &mut cache);

    let mut parts: HashSet<Vec<usize>> = HashSet::new();
    let mut lengths: HashMap<Vec<usize>, f64> = HashMap::new();

    for (node_id, leaves_vec) in cache.iter() {
        if *node_id == root_id { continue; }
        if leaves_vec.len() <= 1 { continue; } // ignore trivial bipartitions
        parts.insert(leaves_vec.clone());
        let node = tree.get(node_id)?;
        // Some BEAST trees may omit branch lengths on certain edges; treat missing as 0.0
        if let Some(len) = node.parent_edge { lengths.insert(leaves_vec.clone(), len); } else { lengths.insert(leaves_vec.clone(), 0.0); }
    }

    let mut root_children: HashSet<Vec<usize>> = HashSet::new();
    let root = tree.get(&root_id)?;
    for &ch in &root.children { if let Some(v) = cache.get(&ch) { root_children.insert(v.clone()); } }

    Ok(PartitionsSnapshot { parts, lengths, root_children, rooted })
}

fn rf_from_snap(a: &PartitionsSnapshot, b: &PartitionsSnapshot) -> usize {
    let inter = a.parts.intersection(&b.parts).count();
    let rf = a.parts.len() + b.parts.len() - 2 * inter;
    let same_root = a.root_children == b.root_children;
    if a.rooted && b.rooted && rf != 0 && !same_root { rf + 2 } else { rf }
}

fn weighted_rf_from_snap(a: &PartitionsSnapshot, b: &PartitionsSnapshot) -> f64 {
    let mut dist = 0.0;
    for p in a.parts.iter() {
        if b.parts.contains(p) {
            let la = a.lengths.get(p).cloned().unwrap_or(0.0);
            let lb = b.lengths.get(p).cloned().unwrap_or(0.0);
            dist += (la - lb).abs();
        } else { dist += a.lengths.get(p).cloned().unwrap_or(0.0); }
    }
    for p in b.parts.iter() { if !a.parts.contains(p) { dist += b.lengths.get(p).cloned().unwrap_or(0.0); } }
    dist
}

fn kf_from_snap(a: &PartitionsSnapshot, b: &PartitionsSnapshot) -> f64 {
    let mut sum = 0.0;
    for p in a.parts.iter() {
        if b.parts.contains(p) {
            let la = a.lengths.get(p).cloned().unwrap_or(0.0);
            let lb = b.lengths.get(p).cloned().unwrap_or(0.0);
            sum += (la - lb).powi(2);
        } else { sum += a.lengths.get(p).cloned().unwrap_or(0.0).powi(2); }
    }
    for p in b.parts.iter() { if !a.parts.contains(p) { sum += b.lengths.get(p).cloned().unwrap_or(0.0).powi(2); } }
    sum.sqrt()
}

// Public custom functions (sequential)
pub fn robinson_foulds_custom(a: &PhyloTree, b: &PhyloTree) -> Result<usize, phylotree::tree::TreeError> { let sa = snapshot_tree_partitions(a)?; let sb = snapshot_tree_partitions(b)?; Ok(rf_from_snap(&sa, &sb)) }
pub fn weighted_robinson_foulds_custom(a: &PhyloTree, b: &PhyloTree) -> Result<f64, phylotree::tree::TreeError> { let sa = snapshot_tree_partitions(a)?; let sb = snapshot_tree_partitions(b)?; Ok(weighted_rf_from_snap(&sa,&sb)) }
pub fn kuhner_felsenstein_custom(a: &PhyloTree, b: &PhyloTree) -> Result<f64, phylotree::tree::TreeError> { let sa = snapshot_tree_partitions(a)?; let sb = snapshot_tree_partitions(b)?; Ok(kf_from_snap(&sa,&sb)) }

// Sequential pairwise helpers (for benchmarking)
pub fn compute_pairwise_rf_custom_seq(trees: &[PhyloTree]) -> Result<Vec<(usize,usize,usize)>, phylotree::tree::TreeError> {
    let snaps: Vec<PartitionsSnapshot> = trees.iter().map(|t| snapshot_tree_partitions(t)).collect::<Result<_,_>>()?;
    let n = snaps.len();
    let mut out = Vec::with_capacity(n.saturating_mul(n.saturating_sub(1))/2);
    for i in 0..n { for j in (i+1)..n { out.push((i,j, rf_from_snap(&snaps[i], &snaps[j]))); } }
    Ok(out)
}
pub fn compute_pairwise_weighted_rf_custom_seq(trees: &[PhyloTree]) -> Result<Vec<(usize,usize,f64)>, phylotree::tree::TreeError> {
    let snaps: Vec<PartitionsSnapshot> = trees.iter().map(|t| snapshot_tree_partitions(t)).collect::<Result<_,_>>()?;
    let n = snaps.len();
    let mut out = Vec::with_capacity(n.saturating_mul(n.saturating_sub(1))/2);
    for i in 0..n { for j in (i+1)..n { out.push((i,j, weighted_rf_from_snap(&snaps[i], &snaps[j]))); } }
    Ok(out)
}
pub fn compute_pairwise_kf_custom_seq(trees: &[PhyloTree]) -> Result<Vec<(usize,usize,f64)>, phylotree::tree::TreeError> {
    let snaps: Vec<PartitionsSnapshot> = trees.iter().map(|t| snapshot_tree_partitions(t)).collect::<Result<_,_>>()?;
    let n = snaps.len();
    let mut out = Vec::with_capacity(n.saturating_mul(n.saturating_sub(1))/2);
    for i in 0..n { for j in (i+1)..n { out.push((i,j, kf_from_snap(&snaps[i], &snaps[j]))); } }
    Ok(out)
}

// Parallel pairwise RF using snapshots (safe parallelization)
pub fn compute_pairwise_rf_custom_parallel(trees: &[PhyloTree]) -> Result<Vec<(usize,usize,usize)>, phylotree::tree::TreeError> {
    let snaps: Vec<PartitionsSnapshot> = trees
        .iter()
        .map(|t| snapshot_tree_partitions(t))
        .collect::<Result<_, _>>()?;
    let n = snaps.len();
    let chunks: Vec<Vec<(usize, usize, usize)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let si = &snaps[i];
            let mut v = Vec::with_capacity(n.saturating_sub(i + 1));
            for j in (i + 1)..n {
                let sj = &snaps[j];
                v.push((i, j, rf_from_snap(si, sj)));
            }
            v
        })
        .collect();
    Ok(chunks.into_iter().flatten().collect())
}

// =============================
// Bitset-based snapshots (Rustacean-friendly: fast, compact, mergeable)
// =============================

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Bitset(pub Vec<u64>);

impl Bitset {
    fn zeros(words: usize) -> Self { Bitset(vec![0u64; words]) }
    #[inline]
    fn set(&mut self, idx: usize) {
        let word = idx >> 6; // /64
        let bit = idx & 63;
        self.0[word] |= 1u64 << bit;
    }
    #[inline]
    fn or_assign(&mut self, other: &Bitset) {
        for (a, b) in self.0.iter_mut().zip(&other.0) { *a |= *b; }
    }
}

#[derive(Debug, Clone)]
pub struct BitsetSnapshot {
    pub parts: Vec<Bitset>,      // sorted lexicographically
    pub lengths: Vec<f64>,       // aligned with parts
    pub root_children: Vec<Bitset>, // sorted per child subtree
    pub words: usize,
    pub rooted: bool,
}

fn snapshot_tree_partitions_bitset(tree: &PhyloTree) -> Result<BitsetSnapshot, phylotree::tree::TreeError> {
    use phylotree::tree::TreeError;
    let rooted = tree.is_rooted()?;
    let leaves = tree.get_leaves();
    if leaves.is_empty() { return Err(TreeError::IsEmpty); }
    let words = (leaves.len() + 63) / 64;

    let mut leaf_pos: HashMap<usize, usize> = HashMap::with_capacity(leaves.len());
    for (i, id) in leaves.iter().enumerate() { leaf_pos.insert(*id, i); }

    // DFS returns Bitset of leaf membership
    fn collect_bitset(
        node_id: usize,
        tree: &PhyloTree,
        leaf_pos: &HashMap<usize, usize>,
        words: usize,
        cache: &mut HashMap<usize, Bitset>,
    ) -> Bitset {
        if let Some(v) = cache.get(&node_id) { return v.clone(); }
        let node = tree.get(&node_id).expect("valid node");
        if node.children.is_empty() {
            let mut bs = Bitset::zeros(words);
            let idx = *leaf_pos.get(&node_id).expect("leaf mapped");
            bs.set(idx);
            cache.insert(node_id, bs.clone());
            return bs;
        }
        let mut agg = Bitset::zeros(words);
        for &ch in &node.children { let sub = collect_bitset(ch, tree, leaf_pos, words, cache); agg.or_assign(&sub); }
        cache.insert(node_id, agg.clone());
        agg
    }

    let mut cache: HashMap<usize, Bitset> = HashMap::new();
    let root_id = tree.get_root()?;
    let _ = collect_bitset(root_id, tree, &leaf_pos, words, &mut cache);

    let mut parts: Vec<Bitset> = Vec::new();
    let mut lengths: Vec<f64> = Vec::new();
    for (node_id, bs) in cache.iter() {
        if *node_id == root_id { continue; }
        // Ignore trivial partitions (<=1 leaf)
        // Fast popcount: sum of popcounts over words
        let mut pop = 0usize;
        for w in &bs.0 { pop += w.count_ones() as usize; if pop > 1 { break; } }
        if pop <= 1 { continue; }
        parts.push(bs.clone());
        let node = tree.get(node_id)?;
        let len = node.parent_edge.unwrap_or(0.0);
        lengths.push(len);
    }
    // Sort by bitset to allow O(m+n) merges
    let mut idx: Vec<usize> = (0..parts.len()).collect();
    idx.sort_unstable_by(|&i, &j| parts[i].cmp(&parts[j]));
    let mut parts_sorted = Vec::with_capacity(parts.len());
    let mut lengths_sorted = Vec::with_capacity(lengths.len());
    for k in idx { parts_sorted.push(parts[k].clone()); lengths_sorted.push(lengths[k]); }

    // Root children bitsets for rooted adjustment
    let root = tree.get(&root_id)?;
    let mut root_children: Vec<Bitset> = Vec::new();
    for &ch in &root.children {
        if let Some(bs) = cache.get(&ch) { root_children.push(bs.clone()); }
    }
    root_children.sort_unstable();

    Ok(BitsetSnapshot { parts: parts_sorted, lengths: lengths_sorted, root_children, words, rooted })
}

// Merge helpers for two sorted snapshots
fn rf_from_bitsets(a: &BitsetSnapshot, b: &BitsetSnapshot) -> usize {
    let mut i = 0usize; let mut j = 0usize; let mut inter = 0usize;
    while i < a.parts.len() && j < b.parts.len() {
        match a.parts[i].cmp(&b.parts[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => { inter += 1; i += 1; j += 1; }
        }
    }
    let rf = a.parts.len() + b.parts.len() - 2 * inter;
    let rooted_adj = a.rooted && b.rooted && rf != 0 && a.root_children != b.root_children;
    if rooted_adj { rf + 2 } else { rf }
}

fn weighted_rf_from_bitsets(a: &BitsetSnapshot, b: &BitsetSnapshot) -> f64 {
    let mut i = 0usize; let mut j = 0usize; let mut dist = 0.0;
    while i < a.parts.len() && j < b.parts.len() {
        match a.parts[i].cmp(&b.parts[j]) {
            std::cmp::Ordering::Less => { dist += a.lengths[i]; i += 1; }
            std::cmp::Ordering::Greater => { dist += b.lengths[j]; j += 1; }
            std::cmp::Ordering::Equal => { dist += (a.lengths[i] - b.lengths[j]).abs(); i += 1; j += 1; }
        }
    }
    while i < a.parts.len() { dist += a.lengths[i]; i += 1; }
    while j < b.parts.len() { dist += b.lengths[j]; j += 1; }
    dist
}

fn kf_from_bitsets(a: &BitsetSnapshot, b: &BitsetSnapshot) -> f64 {
    let mut i = 0usize; let mut j = 0usize; let mut sum = 0.0_f64;
    while i < a.parts.len() && j < b.parts.len() {
        match a.parts[i].cmp(&b.parts[j]) {
            std::cmp::Ordering::Less => { sum += a.lengths[i] * a.lengths[i]; i += 1; }
            std::cmp::Ordering::Greater => { sum += b.lengths[j] * b.lengths[j]; j += 1; }
            std::cmp::Ordering::Equal => { let d = a.lengths[i] - b.lengths[j]; sum += d * d; i += 1; j += 1; }
        }
    }
    while i < a.parts.len() { sum += a.lengths[i] * a.lengths[i]; i += 1; }
    while j < b.parts.len() { sum += b.lengths[j] * b.lengths[j]; j += 1; }
    sum.sqrt()
}

// Public: build bitset snapshots once for reuse
pub fn build_bitset_snapshots(trees: &[PhyloTree]) -> Result<Vec<BitsetSnapshot>, phylotree::tree::TreeError> {
    trees.iter().map(|t| snapshot_tree_partitions_bitset(t)).collect()
}

// Pairwise (bitset) sequential
pub fn compute_pairwise_rf_bitset_seq(snaps: &[BitsetSnapshot]) -> Vec<(usize, usize, usize)> {
    let n = snaps.len();
    let mut out = Vec::with_capacity(n.saturating_mul(n.saturating_sub(1)) / 2);
    for i in 0..n { for j in (i + 1)..n { out.push((i, j, rf_from_bitsets(&snaps[i], &snaps[j]))); } }
    out
}
pub fn compute_pairwise_weighted_rf_bitset_seq(snaps: &[BitsetSnapshot]) -> Vec<(usize, usize, f64)> {
    let n = snaps.len();
    let mut out = Vec::with_capacity(n.saturating_mul(n.saturating_sub(1)) / 2);
    for i in 0..n { for j in (i + 1)..n { out.push((i, j, weighted_rf_from_bitsets(&snaps[i], &snaps[j]))); } }
    out
}
pub fn compute_pairwise_kf_bitset_seq(snaps: &[BitsetSnapshot]) -> Vec<(usize, usize, f64)> {
    let n = snaps.len();
    let mut out = Vec::with_capacity(n.saturating_mul(n.saturating_sub(1)) / 2);
    for i in 0..n { for j in (i + 1)..n { out.push((i, j, kf_from_bitsets(&snaps[i], &snaps[j]))); } }
    out
}

// Pairwise (bitset) parallel
pub fn compute_pairwise_rf_bitset_parallel(snaps: &[BitsetSnapshot]) -> Vec<(usize, usize, usize)> {
    let n = snaps.len();
    let chunks: Vec<Vec<(usize, usize, usize)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let si = &snaps[i];
            let mut v = Vec::with_capacity(n.saturating_sub(i + 1));
            for j in (i + 1)..n { v.push((i, j, rf_from_bitsets(si, &snaps[j]))); }
            v
        })
        .collect();
    chunks.into_iter().flatten().collect()
}
pub fn compute_pairwise_weighted_rf_bitset_parallel(snaps: &[BitsetSnapshot]) -> Vec<(usize, usize, f64)> {
    let n = snaps.len();
    let chunks: Vec<Vec<(usize, usize, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let si = &snaps[i];
            let mut v = Vec::with_capacity(n.saturating_sub(i + 1));
            for j in (i + 1)..n { v.push((i, j, weighted_rf_from_bitsets(si, &snaps[j]))); }
            v
        })
        .collect();
    chunks.into_iter().flatten().collect()
}
pub fn compute_pairwise_kf_bitset_parallel(snaps: &[BitsetSnapshot]) -> Vec<(usize, usize, f64)> {
    let n = snaps.len();
    let chunks: Vec<Vec<(usize, usize, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let si = &snaps[i];
            let mut v = Vec::with_capacity(n.saturating_sub(i + 1));
            for j in (i + 1)..n { v.push((i, j, kf_from_bitsets(si, &snaps[j]))); }
            v
        })
        .collect();
    chunks.into_iter().flatten().collect()
}

// Parallel pairwise Weighted RF using snapshots
pub fn compute_pairwise_weighted_rf_custom_parallel(
    trees: &[PhyloTree],
) -> Result<Vec<(usize, usize, f64)>, phylotree::tree::TreeError> {
    let snaps: Vec<PartitionsSnapshot> = trees
        .iter()
        .map(|t| snapshot_tree_partitions(t))
        .collect::<Result<_, _>>()?;
    let n = snaps.len();
    let chunks: Vec<Vec<(usize, usize, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let si = &snaps[i];
            let mut v = Vec::with_capacity(n.saturating_sub(i + 1));
            for j in (i + 1)..n {
                let sj = &snaps[j];
                v.push((i, j, weighted_rf_from_snap(si, sj)));
            }
            v
        })
        .collect();
    Ok(chunks.into_iter().flatten().collect())
}

// Parallel pairwise Kuhner–Felsenstein (branch score) using snapshots
pub fn compute_pairwise_kf_custom_parallel(
    trees: &[PhyloTree],
) -> Result<Vec<(usize, usize, f64)>, phylotree::tree::TreeError> {
    let snaps: Vec<PartitionsSnapshot> = trees
        .iter()
        .map(|t| snapshot_tree_partitions(t))
        .collect::<Result<_, _>>()?;
    let n = snaps.len();
    let chunks: Vec<Vec<(usize, usize, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let si = &snaps[i];
            let mut v = Vec::with_capacity(n.saturating_sub(i + 1));
            for j in (i + 1)..n {
                let sj = &snaps[j];
                v.push((i, j, kf_from_snap(si, sj)));
            }
            v
        })
        .collect();
    Ok(chunks.into_iter().flatten().collect())
}

#[cfg(test)]
mod tests {}

#[cfg(test)]
mod custom_tests {
    use super::*;
    fn ct1() -> PhyloTree { PhyloTree::from_newick("(A:0.1,(B:0.2,(C:0.3,D:0.4):0.5):0.6,(E:0.7,(F:0.8,G:0.9):0.95):1.0,H:1.1)R;").unwrap() }
    fn ct2() -> PhyloTree { PhyloTree::from_newick("((A:0.1,B:0.2):0.25,(C:0.3,(D:0.4,E:0.5):0.55):0.6,(F:0.8,(G:0.9,H:1.0):0.85):0.9)R;").unwrap() }
    fn ct3() -> PhyloTree { PhyloTree::from_newick("(A:0.1,B:0.2,(C:0.3,(D:0.4,E:0.5):0.55):0.6,(F:0.8,(G:0.9,H:1.0):0.85):0.9)R;").unwrap() }
    fn ct4() -> PhyloTree { PhyloTree::from_newick("((A:0.1,(B:0.2,C:0.3):0.35):0.4,(D:0.4,E:0.5,(F:0.8,(G:0.9,H:1.0):0.85):0.9):0.95)R;").unwrap() }

    #[test]
    fn custom_rf_matches_library() {
        let t1 = ct1(); let t2 = ct2();
        let lib = t1.robinson_foulds(&t2).unwrap();
        let custom = robinson_foulds_custom(&t1,&t2).unwrap();
        assert_eq!(lib, custom);
    }

    #[test]
    fn custom_weighted_matches_library() {
        let t1 = ct3(); let t2 = ct4();
        let lib = t1.weighted_robinson_foulds(&t2).unwrap();
        let custom = weighted_robinson_foulds_custom(&t1,&t2).unwrap();
        assert!((lib - custom).abs() < 1e-9, "weighted mismatch lib={lib} custom={custom}");
    }

    #[test]
    fn custom_kf_matches_library() {
        let t1 = ct1(); let t2 = ct4();
        let lib = t1.khuner_felsenstein(&t2).unwrap();
        let custom = kuhner_felsenstein_custom(&t1,&t2).unwrap();
        assert!((lib - custom).abs() < 1e-9);
    }

    #[test]
    fn parallel_rf_custom() {
        let trees = vec![ct1(), ct2(), ct3(), ct4()];
        let mut seq: Vec<(usize,usize,usize)> = Vec::new();
        for i in 0..trees.len() {
            for j in (i+1)..trees.len() {
                seq.push((i, j, robinson_foulds_custom(&trees[i], &trees[j]).unwrap()));
            }
        }
        let par = compute_pairwise_rf_custom_parallel(&trees).unwrap();
        assert_eq!(seq, par);
    }

    #[test]
    fn parallel_weighted_and_kf_custom() {
        let trees = vec![ct1(), ct2(), ct3(), ct4()];
        // Weighted RF
        let mut seq_w: Vec<(usize,usize,f64)> = Vec::new();
        for i in 0..trees.len() {
            for j in (i+1)..trees.len() {
                seq_w.push((i, j, weighted_robinson_foulds_custom(&trees[i], &trees[j]).unwrap()));
            }
        }
        let par_w = compute_pairwise_weighted_rf_custom_parallel(&trees).unwrap();
        assert_eq!(seq_w.len(), par_w.len());
        for (a,b) in seq_w.iter().zip(par_w.iter()) {
            assert_eq!(a.0, b.0);
            assert_eq!(a.1, b.1);
            assert!((a.2 - b.2).abs() < 1e-12);
        }

        // KF
        let mut seq_k: Vec<(usize,usize,f64)> = Vec::new();
        for i in 0..trees.len() {
            for j in (i+1)..trees.len() {
                seq_k.push((i, j, kuhner_felsenstein_custom(&trees[i], &trees[j]).unwrap()));
            }
        }
        let par_k = compute_pairwise_kf_custom_parallel(&trees).unwrap();
        assert_eq!(seq_k.len(), par_k.len());
        for (a,b) in seq_k.iter().zip(par_k.iter()) {
            assert_eq!(a.0, b.0);
            assert_eq!(a.1, b.1);
            assert!((a.2 - b.2).abs() < 1e-12);
        }
    }
}



