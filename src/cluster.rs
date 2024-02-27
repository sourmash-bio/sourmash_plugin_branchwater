use anyhow::{Context, Result};
use rustworkx_core::connectivity::connected_components;
use rustworkx_core::petgraph::graph::{NodeIndex, UnGraph};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

use crate::utils::MultiSearchResult;

// potential todo:
// - eval DiGraph for directed similarity info (e.g. input containment_A, containment_B independently)
// - explore if collect-first, add edges second style parallelization is worthwhile

fn build_graph(
    file_path: &str,
    similarity_measure: &str,
    similarity_threshold: f64,
) -> Result<(UnGraph<String, f64>, HashMap<String, NodeIndex>)> {
    let mut reader = csv::Reader::from_path(file_path).context("Failed to open CSV file")?;
    let mut name_to_node: HashMap<String, NodeIndex> = HashMap::new();
    let mut graph = UnGraph::<String, f64>::new_undirected();

    for result in reader.deserialize::<MultiSearchResult>() {
        let record = result.map_err(|e| anyhow::anyhow!("Error deserializing record: {}", e))?;

        // ignore self-matches reported via multisearch
        if record.query_name == record.match_name {
            continue;
        }

        let similarity = match similarity_measure {
            "containment" => record.containment,
            "max_containment" => record.max_containment,
            "jaccard" => record.jaccard,
            "average_containment_ani" => match record.average_containment_ani {
                Some(value) => value,
                None => {
                    return Err(anyhow::anyhow!(
                        "average_containment_ani is None. Did you estimate ANI?"
                    ))
                }
            },
            "max_containment_ani" => match record.max_containment_ani {
                Some(value) => value,
                None => {
                    return Err(anyhow::anyhow!(
                        "max_containment_ani is None. Did you estimate ANI?"
                    ))
                }
            },
            _ => {
                return Err(anyhow::anyhow!(
                    "Invalid similarity measure: {}",
                    similarity_measure
                ))
            } // should not happen
        };

        let node1 = *name_to_node
            .entry(record.query_name.clone())
            .or_insert_with(|| graph.add_node(record.query_name.clone()));
        let node2 = *name_to_node
            .entry(record.match_name.clone())
            .or_insert_with(|| graph.add_node(record.match_name.clone()));

        if similarity >= similarity_threshold {
            graph.add_edge(node1, node2, similarity);
        }
    }

    if graph.node_count() == 0 {
        bail!("No nodes added to graph.")
    }

    if graph.edge_count() == 0 {
        bail!("Graph has nodes but no edges were added.");
    }

    Ok((graph, name_to_node))
}

pub fn cluster(
    pairwise_csv: String,
    output_clusters: String,
    similarity_column: String,
    similarity_threshold: f64,
    cluster_sizes: Option<String>,
) -> Result<()> {
    let (graph, name_to_node) =
        match build_graph(&pairwise_csv, &similarity_column, similarity_threshold) {
            Ok(result) => result,
            Err(e) => {
                eprintln!("Error: {:?}", e); // print the underlying error.
                bail!("Failed to build graph.");
            }
        };
    let components = connected_components(&graph);

    // HashMap to count cluster sizes
    let mut size_counts: HashMap<usize, usize> = HashMap::new();

    // Open file for components + names
    let mut file = File::create(output_clusters).context("Failed to create output file")?;

    // write header
    writeln!(file, "cluster,nodes").context("Failed to write header to output file")?;
    // for each component, find corresponding node names + write to file
    for (i, component) in components.iter().enumerate() {
        let component_name = format!("Component_{}", i + 1);
        let node_names: Vec<String> = component
            .iter()
            .filter_map(|node_id| {
                name_to_node.iter().find_map(|(name, &id)| {
                    if id == *node_id {
                        let ident: &str = name.split(' ').next().unwrap();
                        Some(ident.to_owned())
                    } else {
                        None
                    }
                })
            })
            .collect();

        let node_names_str = node_names.join(";");

        writeln!(file, "{},{}", component_name, node_names_str).context(format!(
            "Failed to write component {} to output file",
            i + 1
        ))?;

        // add cluster to aggregated counts
        let count = size_counts.entry(component.len()).or_insert(0);
        *count += 1;
    }

    // write the sizes and counts
    if let Some(sizes_file) = cluster_sizes {
        let mut cluster_size_file =
            File::create(sizes_file).context("Failed to create cluster size file")?;
        writeln!(cluster_size_file, "cluster_size,count")
            .context("Failed to write header to cluster size file")?;
        for (size, count) in size_counts {
            writeln!(cluster_size_file, "{},{}", size, count)
                .context("Failed to write size count to cluster size file")?;
        }
    }

    Ok(())
}
