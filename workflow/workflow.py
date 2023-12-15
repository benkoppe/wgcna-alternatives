import pandas as pd
import numpy as np
from pathlib import Path
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

pandas2ri.activate()

# ------------- SETTINGS --------------
dir = Path(__file__).parent

# Path to gene expression file and its separator
# file must:
#   - have samples as columns and genes as rows
#   - fill first row with sample IDs
#   - fill first column with gene IDs
expression_fpath = dir / "data/expression.txt"
expression_separator = "\t"

# Settings for filtering by sample information
#
# Samples must be listed in the same order as in the expression file,
# and the file must contain a named row: a so-called 'match field'
# that contains some string that must be matched. In this case,
# that string is "areaX". All samples not containing the string in the
# match field will be dropped.
#
# To disable, set sample information fpath to None
sample_information_fpath = dir / "data/sample info.txt"
sample_information_separator = "\t"

match_field_rowname = "Sample_number"
match_value = "areaX"

# Settings for filtering NA values and by sample sum thresholds.
#
# Replaces NA values with 0, and removes samples with a sum outside a certain
# order of magnitude of the median of all other sample sums. The
# order_of_magnitude_threshold setting sets the order needed for a sample to be removed.
#
# To disable, set order_of_magnitude_threshold to None
order_of_magnitude_threshold = 2

# Settings for filtering by gene standard deviation thresholds
#
# All genes with values greater than two standard deviations outside the other data in the sample
# causes the entire gene to be removed. This controls for abnormal and poor data
#
# To disable, set num_deviations_threshold to None
num_deviations_threshold = 2

# Settings for experimental GWENA filtering
#
# Uses built-in filtering in the GWENA package. Currently, this filter is likely not being used properly,
# as nothing is every filtered out
gwena_filter_enable = True

# Save directory
#
# All plots and reports will be saved to this directory
save_dir = dir / "output"


# ------------- BEGIN SCRIPT --------------
def get_df_size(df):
    """Returns the number of columns and rows in a df"""
    return len(df.columns), len(df)


def print_df_changes(step_name, initial_size, final_size):
    """Prints the change in columns and rows between two states of a df"""
    print(f"{step_name}:")
    print(
        f"\tColumns changed from {initial_size[0]} to {final_size[0]}, rows changed from {initial_size[1]} to {final_size[1]}"
    )


def main():
    # -- Loading DataFrame (df) --
    df = pd.read_csv(expression_fpath, sep=expression_separator)
    df.columns = df.columns.str.strip()

    # -- Filtering by sample information --
    if sample_information_fpath is not None:
        # load file
        labels = pd.read_csv(sample_information_fpath, sep=sample_information_separator)
        labels.set_index(labels.columns[0], inplace=True)
        labels = labels.transpose()

        initial_size = get_df_size(df)  # stores the size of df before sample filtering

        # loop for matches and drop columns
        for col_name in df:
            col_type = labels.loc[labels[match_field_rowname] == col_name].index[0]

            if match_value not in col_type:
                df = df.drop(col_name, axis=1)

        print_df_changes(
            "Sample information filtering", initial_size, get_df_size(df)
        )  # prints change in df

    # -- NA/Sum Filtering --
    if order_of_magnitude_threshold is not None:
        initial_size = get_df_size(df)  # stores the size of df before sample filtering

        df.fillna(0)  # replace all NA and NaN values with 0

        # compute sums and medians
        sums = {col_name: np.sum(col_data) for col_name, col_data in df.items()}
        median = np.median(list(sums.values()))

        # loop for sums outside order of magnitude threshold
        for col_name, _ in df.items():
            sum = sums[col_name]
            ratio = sum / median

            if ratio > (10**order_of_magnitude_threshold) or ratio < (
                10**-order_of_magnitude_threshold
            ):
                df = df.drop(col_name, axis=1)

        print_df_changes(
            "NA/Sum filtering", initial_size, get_df_size(df)
        )  # prints change in df

    # -- Standard Deviation Filtering --
    if num_deviations_threshold is not None:
        initial_size = get_df_size(df)  # stores the size of df before sample filtering

        # compute means and standard deviations
        means_and_devs = {
            col_name: (np.mean(col_data), np.std(col_data))
            for col_name, col_data in df.items()
        }

        # loop for values outside standard deviation threshold
        for col_name, col_data in df.items():
            mean, dev = means_and_devs[col_name]

            upper_bound = mean + (dev * num_deviations_threshold)
            lower_bound = mean - (dev * num_deviations_threshold)

            for idx, item in col_data.items():
                if (item < lower_bound) or (item > upper_bound):
                    # mark cells that meet criteria
                    df.at[idx, col_name] = np.nan

        # remove marked rows
        df = df.dropna()  # FIXME: Data hazard when NA filtering is disabled.

        print_df_changes(
            "Standard deviation filtering", initial_size, get_df_size(df)
        )  # prints change in df

    # -- EXPERIMENTAL: GWENA Filtering --
    if gwena_filter_enable:
        initial_size = get_df_size(df)  # stores the size of df before sample filtering

        limma = importr("limma")
        df_norm = limma.normalizeBetweenArrays(df, method="quantile")
        df[:] = df_norm

        gwena = importr("GWENA")  # imports GWENA package
        df_filter = gwena.filter_low_var(df)  # calls GWENA filtering function
        df = pd.DataFrame(
            df_filter
        ).transpose()  # converts result back to pandas DataFrame

        print_df_changes(
            "GWENA filtering", initial_size, get_df_size(df)
        )  # prints change in df

    # -- Running CEMiTool --
    cemitool = importr("CEMiTool")  # import CEMiTool
    doParallel = importr("doParallel")  # import doParallel

    doParallel.registerDoParallel(cores=8)  # configure doParallel
    print("Cemitool and parallel initialized")

    cem = cemitool.cemitool(df, plot=True, verbose=True)  # run CEMiTool

    # -- Save Data --
    cemitool.generate_report(cem, force=True, directory=str(save_dir / "report"))
    cemitool.diagnostic_report(cem, force=True, directory=str(save_dir / "diagnostic"))
    cemitool.write_files(cem, force=True, directory=str(save_dir / "files"))
    cemitool.save_plots(cem, "all", force=True, directory=str(save_dir / "plots"))


if __name__ == "__main__":
    main()
