import argparse
import numpy as np
import pandas as pd
from astropy.table import Table


def get_one_lc_slopes(df, timerange=0, verbose=False, sanity_plots=False):
    """
    
    Args:
        df (pd.DataFrame): containing m, dm
        timerange (int): days to take into account (0=all)
        verbose (Boolean): print verbose
        sanity_plots (Boolean): plot a light-curve and slopes
    Returns:
        dic_out (dict): dictionary with values for the light-curve
    """

    dic_out = {}

    if timerange != 0:
        df_sel = df[df["MJD"] > (df["MJD"].max() - timerange)]
    else:
        df_sel = df

    df_sel_new = df_sel[int(len(df) / 2) :]
    df_sel_old = df_sel[: int(len(df) / 2)]

    dic_out["mad_slope"] = df_sel_new["m"].median() - df_sel_old["m"].median()

    dic_out["weighted_average_slope"] = np.sum(
        df_sel_new["m"] / df_sel_new["dm"]
    ) / np.sum(1 / df_sel_new["dm"]) - np.sum(
        df_sel_old["m"] / df_sel_old["dm"]
    ) / np.sum(
        1 / df_sel_old["dm"]
    )

    dic_out["last_measurement_slope"] = df["m"].iloc[-1] - df["m"].iloc[-2]

    # rate per day
    norm_days = df_sel.MJD.max() - df_sel.MJD.min()
    dic_out["mad_slope_rate"] = dic_out["mad_slope"] / norm_days
    dic_out["weighted_average_slope_rate"] = (
        dic_out["weighted_average_slope"] / norm_days
    )

    if verbose == True:
        print(dic_out)
    if sanity_plots == True:
        import matplotlib.pyplot as plt

        plt.errorbar(df.MJD, df.m, yerr=df.dm, fmt="o")
        days = np.arange(df.MJD.min(), df.MJD.max() + 0.1, 0.1)
        mad_line = df.m[0] + dic_out["mad_slope"] * (days - df.MJD.min())
        plt.plot(days, mad_line, label="mad")
        weighted_average_line = df.m[0] + dic_out["weighted_average_slope"] * (
            days - df.MJD.min()
        )
        plt.plot(days, weighted_average_line, label="wa")

        mad_line = df.m[0] + dic_out["mad_slope_rate"] * (days - df.MJD.min())
        plt.plot(days, mad_line, label="mad rate")
        weighted_average_line = df.m[0] + dic_out["weighted_average_slope_rate"] * (
            days - df.MJD.min()
        )
        plt.plot(days, weighted_average_line, label="wa rate")
        plt.legend()
        plt.gca().invert_yaxis()
        plt.savefig("tmp.png")

    return dic_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process candidates metadata")

    parser.add_argument(
        "--fname",
        default="./test_files/lc.txt",
        help="Candidate metadata (coords) from Mary",
    )
    parser.add_argument(
        "--path_out", default="./dump", help="Path to output enriched metadata",
    )
    parser.add_argument(
        "--outname",
        default="transients_coo_enriched",
        help="Output name enriched metadata (no extension)",
    )

    parser.add_argument(
        "--test", action="store_true", help="one file processed only",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.fname, index_col=False, delimiter=" ")

    dic_out = get_one_lc_slopes(df, verbose=True)

