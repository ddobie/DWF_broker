""" DWF brokerage processing
- Ingest candidates coordinates & cross-match to catalogues
- ingest photometry & compute fatures
"""
import os
import argparse
import pandas as pd
from utils import xmatch
from utils import data_utils as du

# TO DO
# argparse
# pickle out?
# USNO!

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process candidates metadata")

    parser.add_argument(
        "--fname",
        default="./test_files/transients_coo.txt",
        help="Candidate metadata (coords) from Mary",
    )
    parser.add_argument(
        "--path_out", default="./dump", help="Path to output enriched metadata",
    )
    parser.add_argument(
        "--outname",
        default="./transients_coo_enriched",
        help="Output name enriched metadata (no extension)",
    )

    parser.add_argument(
        "--test", action="store_true", help="one file processed only",
    )
    args = parser.parse_args()

    df = du.read_mary_masterlist(args.fname)
    if args.test:
        df = df[:10]
    print(f"Processing {args.fname}: {len(df)} candidates")

    # SIMBAD
    print("SIMBAD xmatch")
    z, sptype, typ, ctlg = xmatch.cross_match_simbad(
        df["index"].to_list(), df["ra"].to_list(), df["dec"].to_list(),
    )
    # save in df
    df["simbad_type"] = typ
    df["simbad_ctlg"] = ctlg
    df["simbad_sptype"] = sptype
    df["simbad_redshift"] = z

    df_gaia = xmatch.cross_match_alerts_raw_generic(
        df["index"].to_list(),
        df["ra"].to_list(),
        df["dec"].to_list(),
        ctlg="vizier:I/355/gaiadr3",
        distmaxarcsec=2,
        columns_to_keep=[
            "DR3Name",
            "RAdeg",
            "DEdeg",
            "Plx",
            "e_Plx",
            "Gmag",
            "e_Gmag",
            "angDist",
        ],
    )
    df_gaia = df_gaia.rename(
        columns={
            col: col + "_gaiaDR3"
            for col in df_gaia.columns
            if col not in ["objectId", "ra", "dec"]
        }
    )
    df_gaia["index"] = df_gaia["objectId"]
    df_gaia["index"] = df_gaia["index"].astype(int)

    df = pd.merge(df, df_gaia, on=["index", "ra", "dec"], how="left")

    print("USNO-A.20 xmatch")
    (source_usno, angdist_usno,) = xmatch.cross_match_usno(
        df["index"].to_list(),
        df["ra"].to_list(),
        df["dec"].to_list(),
        ctlg="vizier:I/252/out",
    )
    df["USNO_source"] = source_usno
    df["USNO_angdist"] = angdist_usno

    # save cross-matched alerts
    os.makedirs(args.path_out, exist_ok=True)
    df.to_csv(f"{args.path_out}/{args.outname}.csv")
    print(
        f"Saved enriched metadata of {len(df)} candidates in {args.path_out}/{args.outname}.csv"
    )
