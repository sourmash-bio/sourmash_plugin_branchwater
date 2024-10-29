import csv


def pretty_print_manysearch(manysearch_csv):
    "Pretty-print the manysearch output."
    with open(manysearch_csv, newline="") as fp:
        r = csv.DictReader(fp)
        rows = list(r)

    rows.sort(key=lambda row: row["query_name"])  # sort by metagenome, for now

    first = True
    for row in rows:
        has_abundance = "average_abund" in row

        #
        # display!
        #

        # displaying first result?
        if first:
            print("")
            print("query             p_genome avg_abund   p_metag   metagenome name")
            print("--------          -------- ---------   -------   ---------------")
            first = False

        f_genome_found = float(row["containment"])
        pct_genome = f"{f_genome_found*100:.1f}"

        if has_abundance:
            n_weighted_found = int(row["n_weighted_found"])
            total_weighted_hashes = int(row["total_weighted_hashes"])
            f_metag_weighted = (
                n_weighted_found / total_weighted_hashes
            )  # results_d['f_match_weighted']
            pct_metag = f"{f_metag_weighted*100:.1f}%"

            avg_abund = float(row["average_abund"])
            avg_abund = f"{avg_abund:.1f}"
        else:
            avg_abund = "N/A"
            pct_metag = "N/A"

        query_name = row["query_name"][:17]
        metag_name = row["match_name"][:17]
        print(
            f"{query_name:<17} {pct_genome:>6}%  {avg_abund:>6}     {pct_metag:>6}     {metag_name}"
        )
