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

    # GAIA DR2
    print("GAIA DR2 xmatch")
    gaia_source, plx, plxerr, gmag, angdist = xmatch.cross_match_gaia(
        df["index"].to_list(),
        df["ra"].to_list(),
        df["dec"].to_list(),
        ctlg="vizier:I/345/gaia2",
    )
    df["gaia_DR2_source"] = gaia_source
    df["gaia_DR2_parallax"] = plx
    df["gaia_DR2_parallaxerr"] = plxerr
    df["gaia_DR2_gmag"] = gmag
    df["gaia_DR2_angdist"] = angdist

    # GAIA eDR3
    print("GAIA eDR3 xmatch")
    gaia_source, plx, plxerr, gmag, angdist = xmatch.cross_match_gaia(
        df["index"].to_list(),
        df["ra"].to_list(),
        df["dec"].to_list(),
        ctlg="vizier:I/350/gaiaedr3",
    )
    df["gaia_eDR3_source"] = gaia_source
    df["gaia_eDR3_parallax"] = plx
    df["gaia_eDR3_parallaxerr"] = plxerr
    df["gaia_eDR3_gmag"] = gmag
    df["gaia_eDR3_angdist"] = angdist

    # save cross-matched alerts
    os.makedirs(args.path_out, exist_ok=True)
    df.to_csv(f"{args.path_out}/{args.outname}.csv")
    print(
        f"Saved enriched metadata of {len(df)} candidates in {args.path_out}/{args.outname}.csv"
    )
