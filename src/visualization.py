import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def plot_variant_histogram(variants, column, bins=10, figsize=(10, 6)):
    """
    Plot a histogram for a specified column of the annotated variants.

    Args:
        variants (pandas.DataFrame): DataFrame containing the annotated variants.
        column (str): Column name to plot the histogram for.
        bins (int, optional): Number of bins for the histogram. Default is 10.
        figsize (tuple, optional): Figure size for the plot. Default is (10, 6).
    """
    plt.figure(figsize=figsize)
    plt.hist(variants[column], bins=bins, edgecolor='black')
    plt.title(f'Histogram of {column}', fontsize=16)
    plt.xlabel(column, fontsize=14)
    plt.ylabel('Count', fontsize=14)
    plt.tight_layout()
    plt.show()

def plot_variant_scatter(variants, x_column, y_column, hue_column=None, figsize=(10, 6)):
    """
    Plot a scatter plot for the annotated variants.

    Args:
        variants (pandas.DataFrame): DataFrame containing the annotated variants.
        x_column (str): Column name for the x-axis.
        y_column (str): Column name for the y-axis.
        hue_column (str, optional): Column name to use for coloring the points. Default is None.
        figsize (tuple, optional): Figure size for the plot. Default is (10, 6).
    """
    plt.figure(figsize=figsize)
    if hue_column:
        sns.scatterplot(data=variants, x=x_column, y=y_column, hue=hue_column)
    else:
        sns.scatterplot(data=variants, x=x_column, y=y_column)
    plt.title(f'{y_column} vs {x_column}', fontsize=16)
    plt.xlabel(x_column, fontsize=14)
    plt.ylabel(y_column, fontsize=14)
    plt.tight_layout()
    plt.show()