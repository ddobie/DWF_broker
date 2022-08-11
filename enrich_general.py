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

    # df = du.read_mary_masterlist(args.fname)
    df = pd.read_csv(args.fname, delimiter="\t")

    def convert_to_degrees(row, key="RA"):
        from astropy import units as u
        from astropy.coordinates import SkyCoord

        c = SkyCoord(row.RA, row.DEC, unit=(u.hourangle, u.deg))
        if key == "RA":
            return c.ra.degree
        elif key == "DEC":
            return c.dec.degree

    df["ra"] = df.apply(
        lambda row: row.RA
        if ":" not in str(row.RA)
        else convert_to_degrees(row, key="RA"),
        axis=1,
    )
    df["dec"] = df.apply(
        lambda row: row.DEC
        if ":" not in str(row.DEC)
        else convert_to_degrees(row, key="DEC"),
        axis=1,
    )
    df["index"] = df.index.values

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
