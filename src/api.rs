// //! Python binding layer (optional). Provides minimal functions for tree loading.
// use pyo3::prelude::*;
// use crate::{read_beast_trees, compute_pairwise_robinson_foulds};

// #[pyfunction]
// #[pyo3(signature = (path, burnin_trees=0, burnin_states=0, use_real_taxa=true))]
// fn py_read_beast_trees(
//     path: String,
//     burnin_trees: usize,
//     burnin_states: i64,
//     use_real_taxa: bool,
// ) -> PyResult<Vec<String>> {
//     let trees = read_beast_trees(&path, burnin_trees, burnin_states, use_real_taxa);
//     Ok(trees.iter().map(|t| format!("{t}")).collect())
// }

// #[pyfunction]
// #[pyo3(signature = (path, burnin_trees=0, burnin_states=0, use_real_taxa=true))]
// fn py_pairwise_rf(
//     path: String,
//     burnin_trees: usize,
//     burnin_states: i64,
//     use_real_taxa: bool,
// ) -> PyResult<Vec<(usize, usize, f64)>> {
//     let trees = read_beast_trees(&path, burnin_trees, burnin_states, use_real_taxa);
//     Ok(compute_pairwise_robinson_foulds(&trees))
// }

// #[pymodule]
// fn rust_python_tree_distances_api(_py: Python, m: &PyModule) -> PyResult<()> {
//     m.add_function(wrap_pyfunction!(py_read_beast_trees, m)?)?;
//     m.add_function(wrap_pyfunction!(py_pairwise_rf, m)?)?;
//     Ok(())
// }
