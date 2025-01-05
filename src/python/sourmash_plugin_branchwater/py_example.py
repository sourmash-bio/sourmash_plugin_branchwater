from . import sourmash_plugin_branchwater as api

def do_fastmultigather(query_paths,
                       against_paths,
                       threshold_bp,
                       ksize,
                       scaled,
                       moltype,
                       output,
                       save_matches,
                       create_empty_results):
    query_coll = api.api_load_collection(query_paths, ksize, scaled, moltype)
    against_coll = api.api_load_collection(against_paths, ksize, scaled, moltype)

    threshold_hashes = int(threshold_bp / scaled)

    res = against_coll.fastmultigather_against(query_coll,
                                               threshold_hashes,
                                               scaled,
                                               output)

    n_processed, n_skipped, n_failed = res
    print(f"DONE. Processed {n_processed} queries total.")

    if n_skipped:
        print(f"WARNING: skipped {n_skipped} query paths - no compatible signatures.",
              file=sys.stderr)

    if n_failed:
        print(f"WARNING: {n_failed} query paths failed to load. See error messages above.",
              file=sys.stderr)

    return 0
