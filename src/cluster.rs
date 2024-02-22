use anyhow::{Context, Result};
use rustworkx_core::connectivity::connected_components;
use rustworkx_core::petgraph::graph::{NodeIndex, UnGraph};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

use crate::utils::MultiSearchResult;

// todo: eval DiGraph for directed similarity info (e.g. input containment_A, containment_B independently)

fn build_graph(
    file_path: &str,
    similarity_measure: &str,
    similarity_threshold: f64,
) -> Result<(UnGraph<String, f64>, HashMap<String, NodeIndex>)> {
    let mut reader = csv::Reader::from_path(file_path).context("Failed to open CSV file")?;
    let mut graph = UnGraph::<String, f64>::new_undirected();
    let mut name_to_node: HashMap<String, NodeIndex> = HashMap::new();

    for result in reader.deserialize() {
        let record: MultiSearchResult = result.context("Error deserializing CSV record")?;
        let similarity = match similarity_measure {
            // start with containment/jaccard for now.
            "containment" => record.containment,
            "jaccard" => record.jaccard,
            _ => return Err(anyhow::anyhow!("Invalid similarity measure")),
        };

        if similarity >= similarity_threshold {
            let node1 = *name_to_node
                .entry(record.query_name.clone())
                .or_insert_with(|| graph.add_node(record.query_name.clone()));
            let node2 = *name_to_node
                .entry(record.match_name.clone())
                .or_insert_with(|| graph.add_node(record.match_name.clone()));

            graph.add_edge(node1, node2, similarity);
        }
    }

    Ok((graph, name_to_node))
}

pub fn cluster(
    pairwise_csv: String,
    output_clusters: String,
    similarity_column: String,
    similarity_threshold: f64,
) -> Result<()> {
    let (graph, name_to_node) =
        build_graph(&pairwise_csv, &similarity_column, similarity_threshold)
            .context("Failed to build graph")?;

    // Assuming connected_components is defined as shown and returns Vec<HashSet<G::NodeId>>
    let components = connected_components(&graph);

    // Prepare to write the components to a file
    let mut file = File::create(&output_clusters).context("Failed to create output file")?;

    // Now, for each component, we find the corresponding node names and write them to the file
    for (i, component) in components.iter().enumerate() {
        let component_name = format!("Component_{}", i + 1);
        let node_names: Vec<String> = component
            .iter()
            .filter_map(|node_id| {
                // Here we reverse lookup our node_id back to the node name
                name_to_node.iter().find_map(|(name, &id)| {
                    if id == *node_id {
                        Some(name.clone())
                    } else {
                        None
                    }
                })
            })
            .collect();

        writeln!(file, "{}: {:?}", component_name, node_names).context(format!(
            "Failed to write component {} to output file",
            i + 1
        ))?;
    }

    Ok(())
}

// #[cfg(test)]
// mod tests {
//     use crate::pairwise;

//     use super::*;
//     use std::fs::File;
//     use std::io::{BufRead, BufReader, Write};
//     use tempfile::tempdir;

//     #[test]
//     fn test_cluster_jaccard() -> Result<()> {
//         let mut pairwise_csv = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
//         pairwise_csv.push = "python/tests/test-data/cluster.pairwise.csv";

//         let temp_dir = tempdir().unwrap();
//         let output_clusters = temp_dir.path().join("output_clusters.txt");
//         let output_clusters_str = output_clusters.to_str().unwrap().to_string();

//         let similarity_column = String::from("jaccard");
//         let similarity_threshold = 0.5;

//         // Call your cluster function
//         cluster(
//             pairwise_csv,
//             similarity_column,
//             output_clusters_str,
//             similarity_threshold,
//         )?;

//         // Read and verify the output
//         let file = File::open(output_clusters)?;
//         let reader = BufReader::new(file);
//         let mut cluster_count = 0;

//         for line in reader.lines() {
//             let line = line?;
//             if line.starts_with("Component") {
//                 cluster_count += 1;
//             }
//         }

//         // Verify the number of clusters is as expected
//         assert_eq!(cluster_count, 2, "Expected 2 clusters");

//         Ok(())
//     }
// }
