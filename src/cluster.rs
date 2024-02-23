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
            "containment" => record.containment,
            "max_containment" => record.max_containment,
            "jaccard" => record.jaccard,
            "ani" => record.average_containment_ani,
            "max_ani" => record.max_containment_ani,
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

    let components = connected_components(&graph);

    // Prepare to write the components to a file
    let mut file = File::create(&output_clusters).context("Failed to create output file")?;

    // for each component, find corresponding node names + write to file
    for (i, component) in components.iter().enumerate() {
        let component_name = format!("Component_{}", i + 1);
        let node_names: Vec<String> = component
            .iter()
            .filter_map(|node_id| {
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
