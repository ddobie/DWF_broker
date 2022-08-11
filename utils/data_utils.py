""" Utils to process data
"""
import pandas as pd
from astropy.table import Table


def read_mary_masterlist(fname):
    """ Read Mary aggregated coordinate file
    Args:
        fname (string): path + name coordinate file
    Returns:
        df (DataFrame): reformatted file
    """
    dat = Table.read(fname, format="ascii")
    df_tmp = dat.to_pandas()
    # reformatting
    df_tmp = df_tmp.rename(
        columns={
            "col1": "objectId",
            "col2": "ra",
            "col3": "dec",
            "col4": "field",
            "col5": "date",
            "col6": "maryseed",
            "col7": "CCD",
            "col8": "ID_CCD",
        }
    )
    df = df_tmp[1:]

    return df
