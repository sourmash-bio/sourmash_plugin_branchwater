import sys
from . import sourmash_plugin_branchwater as api

def do_fastmultigather(query_path,
                       against_path,
                       threshold_bp,
                       ksize,
                       scaled,
                       moltype,
                       output,
                       save_matches,
                       create_empty_results):
    if api.is_revindex_database(against_path):
        is_rocksdb = True
        against_db = api.BranchRevIndex(against_path)
        # @CTB check scaled, ksize, moltype here? or in a 'select'?
        _, against_scaled = against_db.min_max_scaled()
    else:
        is_rocksdb = False
        against_coll = api.api_load_collection(against_path, ksize, scaled, moltype)
        against_scaled = against_coll.max_scaled()

    try:
        query_coll = api.api_load_collection(query_path, ksize, scaled, moltype)
    except:
        return 1

    if scaled is None:
        query_scaled = query_coll.max_scaled()
        scaled = max(query_scaled, against_scaled)
#    else:
#        if scaled > max(query_scaled, against_scaled):
#            print('ERROR, incompatible scaled')
#            return 1

    threshold_hashes = int(threshold_bp / scaled)

    if is_rocksdb:
        res = against_db.fastmultigather_against(query_coll,
                                                 against_db.selection(), # @CTB
                                                 threshold_hashes,
                                                 output)
    else:
        res = against_coll.fastmultigather_against(query_coll,
                                                   threshold_hashes,
                                                   scaled,
                                                   output,
                                                   save_matches=save_matches)

    n_processed, n_skipped, n_failed = res
    print(f"DONE. Processed {n_processed} queries total.")

    if n_skipped:
        print(f"WARNING: skipped {n_skipped} query paths - no compatible signatures.",
              file=sys.stderr)

    if n_failed:
        print(f"WARNING: {n_failed} query paths failed to load. See error messages above.",
              file=sys.stderr)

    return 0
