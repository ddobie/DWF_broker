import argparse
import pandas as pd

list_simbad_galaxies = [
    "galaxy",
    "Galaxy",
    "EmG",
    "EmissionG",
    "Seyfert",
    "Seyfert_1",
    "Seyfert_2",
    "Seyfert1",
    "Seyfert2",
    "BlueCompG",
    "BlueCompactG",
    "StarburstG",
    "LSB_G",
    "LowSurfBrghtG",
    "HII_G",
    "HIIG",
    "High_z_G",
    "GinPair",
    "GinGroup",
    "GtowardsGroup",
    "BClG",
    "BrightestCG",
    "GinCl",
    "GtowardsCl",
    "PartofG",
    "Compact_Gr_G",
    "InteractingG",
    "PairG",
    "GroupG",
    "ClG",
    "SuperClG",
    "Void",
]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process candidates metadata")

    parser.add_argument(
        "--path_out", default="./dump", help="Path to output enriched metadata",
    )
    parser.add_argument(
        "--fname",
        default="dump/transients_coo_enriched.csv",
        help="Enriched metadata used for filtering",
    )
    parser.add_argument(
        "--nameout", default="transients_coo_filtered.csv", help="Filtered metadata",
    )
    parser.add_argument(
        "--test", action="store_true", help="one file processed only",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.fname)
    if args.test:
        df = df[:10]
    print(f"Processing {len(df)} candidates from {args.fname}")

    #
    # LOGIC
    #

    # Cross-match with SIMBAD is Unknown, Fail or a galaxy
    # eliminates known vairable stars or AGNs, etc
    keepcds = ["Unknown", "Fail"] + list_simbad_galaxies
    cut_simbad = df.simbad_type.isin(keepcds)

    # Cross-match with USNO is Unknown
    # eliminates bright stars in the catalogue
    cut_usno = df["USNO_source"].isin(["Unknown", "Fail"])

    # Cross-match with Gaia is a star
    # Criteria of parallax SNR is needed
    # ! logic is different in this cut as Gaia contains extragalactic objects
    gaia_matched_source = ~df["GaiaDR3_DR3Name"].isin(["Unknown", "Fail"])
    df_gaia_matched_source = df[gaia_matched_source]
    df_stars_wgood_plx = df_gaia_matched_source[
        df_gaia_matched_source["GaiaDR3_Plx"].astype(float)
        / df_gaia_matched_source["GaiaDR3_e_Plx"].astype(float)
        > 5
    ]
    cut_gaia = ~df["objectId"].isin(df_stars_wgood_plx["objectId"].values)

    # "Basic cut"
    df_sel = df[cut_simbad & cut_usno & cut_gaia]

    fname_out = f"{args.path_out}{args.nameout}"
    df_sel.to_csv(fname_out)
    print(f"Saved filtered {len(df_sel)} candidates in {fname_out}")

