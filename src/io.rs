use std::fs;
use std::path::Path;
use phylotree::tree::Tree;
use std::collections::HashMap;

pub fn read_beast_trees<P: AsRef<Path>>(
    path: P,
    burnin_trees: usize,
    burnin_states: usize,
    use_real_taxa: bool,
) -> Vec<Tree> {
    let content = match fs::read_to_string(path.as_ref()) {
        Ok(s) => s,
        Err(e) => { eprintln!("Failed to read {:?}: {e}", path.as_ref()); return Vec::new(); }
    };

    let base_name = path.as_ref()
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown");

    let translate = if use_real_taxa { parse_taxon_block(&content) } else { Default::default() };

    let trees = collect_tree_blocks(&content)
        .into_iter()
        .enumerate()

        //generate tree name & extract state number
        .map(|(idx,tree)| {
            let state = extract_state(&tree.header);
            (idx, tree, state, format!("{base_name}_tree_STATE{state}"))
        })

        // Filter out burn-in trees based on count and/or state number if 0 we don't filter
        .filter(|(idx, _tree, state, _name)| {
                (burnin_trees == 0 && burnin_states == 0) ||
                (burnin_trees > 0 && *idx > burnin_trees) ||
                (burnin_states > 0 && *state > burnin_states)
        })

        // read in the files
        .filter_map(|(idx,tree, _state, _name)| {
            let newick = &tree.body;
            let mut phylo_tree = match phylotree::tree::Tree::from_newick(newick) {
                Ok(t) => t,
                Err(e) => {
                    eprintln!("Failed to parse tree {} at index {}: {}", path.as_ref().display(), idx, e);
                    return None;
                }
            };

            // Rename the leaves with the map
            if use_real_taxa {
                rename_leaf_nodes(&mut phylo_tree, &translate);
            }

            Some(phylo_tree)
        })
        .collect::<Vec<_>>();
    trees
}

fn extract_state(header: &str) -> usize {
    if let Some(start) = header.to_ascii_uppercase().find("STATE_") {
        let num_start = start + 6; // length of "STATE_"
        let rest = &header[num_start..];
        let state = rest.chars()
            .take_while(|c| c.is_ascii_digit())
            .collect::<String>();
        if let Ok(num) = state.parse::<usize>() {
            return num;
        }
    }
    0
}

struct TreeBlock<'a> { header: &'a str, body: String }

fn collect_tree_blocks(content: &str) -> Vec<TreeBlock<'_>> {
    content
        .lines()
        .skip_while(|line| !line.to_ascii_uppercase().starts_with("TREE "))
        .take_while(|line| !line.trim().to_ascii_uppercase().starts_with("END;"))
        .filter_map(|line| {
            let mut parts = line.splitn(2, " = ");
            let header = parts.next()?.trim();
            let body = parts.next()?.trim().to_string();
            Some(TreeBlock { header, body })
        })
        .collect()
}

fn parse_taxon_block(content: &str) -> HashMap<String, String> {
    content
        .lines()
        .skip_while(|line| !line.trim().to_ascii_uppercase().starts_with("TRANSLATE"))
        .skip(1)
        .take_while(|line| !line.trim().to_ascii_uppercase().starts_with(";"))
        // STRUCTURE:
        // 1 '1959.M.CD.59.ZR59',
		// 2 '1960.DRC60A',
        .filter_map(|line| {
            let line = line.trim().trim_end_matches(',');
            let mut parts = line.split_whitespace();
            let id = parts.next()?.to_string();
            let label = parts.next()?.trim_matches('\'').to_string();
            Some((id, label))
        })
        .collect::<HashMap<_, _>>()
}

pub fn rename_leaf_nodes(phylo_tree: &mut Tree, translate: &std::collections::HashMap<String, String>) {
    for leaf_id in phylo_tree.get_leaves() {
        if let Ok(node) = phylo_tree.get_mut(&leaf_id) {
            node.name = node
                .name
                .as_ref()
                .and_then(|n| translate.get(n).cloned());
        }
    }
}