import ternary
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.ndimage import gaussian_filter
import argparse

"""

This program generates a ternary plot of METs (Metabolic Equivalent of Task) values
 for three activity categories: work, travel, and recreation. 

The input data must be in a data_XXX.csv format (where XXX is any 3-letter code) 
and should have the following structure:

5 columns x any number of rows
    ID: Respondent ID
    PtotalMETs: Total METs
    percentage_work: Percentage of METs related to work (in %)
    percentage_travel: Percentage of METs related to travel (in %)
    percentage_recreation: Percentage of METs related to recreation (in %)
    weight: Weight balance coefficient (0-1) for weighting each sample (optional)

The program accepts the following command line options:

    --country (Required): The 3-letter code used in the input file name (data_XXX.csv)
    --filter (Optional): Apply a Gaussian filter to the data
    --mets (Optional): Display the average METs of the data points in each histogram bin as circles
    --average (Optional): Display the overall average value on the ternary plot

"""

NUMBER_OF_NBINS = 11
GAUSSIAN_FILTER_SIGMA = 2
BUBBLE_SIZE_RATIO = 10

def custom_formatter(x, pos):
    return f"{100 * x:.1f}"

def calculate_scatter_radius(df, binx):
    """
    Given a dataframe and the number of histogram bins, this function calculates the scatterplot bubble
    radii and returns a dictionary with the calculated radii.

    Parameters:
    -----------
    df : pandas.DataFrame
        The input dataframe containing the data to be used for calculating the radii of the scatterplot bubbles.
    binx : int
        The number of histogram bins to be used in the calculation.

    Returns:
    --------
    R_dict : dict
        A dictionary with keys representing the unique data points, and values representing the corresponding
        calculated radii for the scatterplot bubbles.
    """

    R_dict = dict()
    for i, x in enumerate(binx):
        for j, y in enumerate(binx):
            for k, z in enumerate(binx):
                if i+j+k == 10:
                    tmp_mean = df[(df['on_x1bin'] == i) & (df['on_x2bin'] == j) & (df['on_x3bin'] == k)]['Total'].mean()

                    df.loc[(df['on_x1bin'] == i) & (df['on_x2bin'] == j) & (df['on_x3bin'] == k), 'wb_total'] = df.loc[(df['on_x1bin'] == i) & (df['on_x2bin'] == j) & (df['on_x3bin'] == k), 'Total'] * df.loc[(df['on_x1bin'] == i) & (df['on_x2bin'] == j) & (df['on_x3bin'] == k), 'weight']
                    tmp_people = df.loc[(df['on_x1bin'] == i) & (df['on_x2bin'] == j) & (df['on_x3bin'] == k), 'weight'].sum()
                    if tmp_people > 0:
                        tmp_wb_mean = df.loc[(df['on_x1bin'] == i) & (df['on_x2bin'] == j) & (df['on_x3bin'] == k), 'wb_total'].sum() / tmp_people
                    else:
                        tmp_wb_mean = 0

                    R_dict[(i,j,k)] = tmp_wb_mean

    return R_dict

def export_graph(df, interp_dict, bubble_dict):

    """
    Exports a ternary plot using the input dataframe, filtered dictionary, bubble size dictionary,
    output filename, and option string.

    Parameters:
    -----------
    df : pandas.DataFrame
        The input dataframe containing the data to be plotted on the ternary graph.
    interp_dict : dict
        The filtered dictionary containing the data after applying a Gaussian filter, or any other
        suitable filter, to smooth out the values.
    bubble_dict : dict
        A dictionary with keys representing the unique data points, and values representing the corresponding
        calculated radii for the scatterplot bubbles.
    filename : str
        The name of the output file to save the generated ternary plot.
    optionstr : str
        A string containing options to customize the appearance and behavior of the ternary plot.

    Returns:
    --------
    None. The function saves the ternary plot to a file with the specified filename.
    """

    fig, tax = ternary.figure(scale=NUMBER_OF_NBINS-1)
    cb_kwargs = {"format": ticker.FuncFormatter(custom_formatter), "shrink": 0.5}
    tax.heatmap(interp_dict, vmin = 0, vmax = 0.01, cbarlabel="% activities", cb_kwargs=cb_kwargs)

    tax.left_axis_label("% activity from recreation", fontsize=10, offset=0.15)
    tax.right_axis_label("% activity from travel", fontsize=10, offset=0.15)
    tax.bottom_axis_label("% activity from work", fontsize=10, offset=0.12)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="white", multiple=1)
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    tax.set_axis_limits({'b': [0, 100], 'l': [0, 100], 'r': [0, 100]})
    tax.get_ticks_from_axis_limits()
    tax.set_custom_ticks(fontsize=10, offset=0.02)

    if args.average:

        ave_work = df['x1'].mean()*10
        ave_travel = df['x2'].mean()*10
        ave_recreation = df['x3'].mean()*10

        p1 = (ave_work, ave_travel, ave_recreation)
        p2 = (ave_work, ave_travel, ave_recreation)
        tax.line(p1, p2, linewidth=3., marker='*', color='gray')
        # print("[nwb]",p1)

        tmp_wb_people = df['weight'].sum()
        wb_ave_work       = 10 * df['x1_wb'].sum() / tmp_wb_people
        wb_ave_travel     = 10 * df['x2_wb'].sum() / tmp_wb_people
        wb_ave_recreation = 10 * df['x3_wb'].sum() / tmp_wb_people

        p1 = (wb_ave_work, wb_ave_travel, wb_ave_recreation)
        p2 = (wb_ave_work, wb_ave_travel, wb_ave_recreation)
        tax.line(p1, p2, linewidth=3., marker='o', color='yellow')
        # print("[wb]",p1)

    if args.mets:
        tax.scatter(bubble_dict, marker='o', c=[b/6000 for a,b in bubble_dict.items()], cmap="gray", alpha=0.5, s=[b/BUBBLE_SIZE_RATIO for a,b in bubble_dict.items()], zorder = 2, edgecolor='black', linewidths=[1]) #, label="Green Diamonds # , color='green'

    tax.scatter([(-3,8.45,  4)], marker='o', facecolors='none', cmap="gray", alpha=0.7, c=[ 600/6000], s=600/BUBBLE_SIZE_RATIO, label="600 METs", edgecolor='black', linewidths=[1])
    tax.scatter([(-3,7.8,4)], marker='o', facecolors='none', cmap="gray", alpha=0.7, c=[1500/6000], s=1500/BUBBLE_SIZE_RATIO, label="1500 METs", edgecolor='black', linewidths=[1])
    tax.scatter([(-3,7,4)], marker='o', facecolors='none', cmap="gray", alpha=0.7, c=[3000/6000], s=3000/BUBBLE_SIZE_RATIO, label="3000 METs", edgecolor='black', linewidths=[1])
    tax.scatter([(-3,6,  4)], marker='o', facecolors='none', cmap="gray", alpha=0.7, c=[6000/6000], s=6000/BUBBLE_SIZE_RATIO, label="6000 METs", edgecolor='black', linewidths=[1])

    fig.text(0.15,0.69,"METs", size='x-small', color='black' ,alpha=1, horizontalalignment='center', fontweight=1000)

    fig.text(0.2225,0.746, "600", color='white' ,alpha=0.8, horizontalalignment='center', fontweight=1000, fontsize=3)
    fig.text(0.204,0.703,"1500", color='white',alpha=0.8, horizontalalignment='center', fontweight=1000, fontsize=4)
    fig.text(0.182,0.648,"3000", color='black',alpha=0.8, horizontalalignment='center', fontweight=1000, fontsize=6)
    fig.text(0.154,0.58,"6000", color='black',alpha=0.8, horizontalalignment='center', fontweight=1000, fontsize=7)

    tax.savefig(filename+"_fig"+optionstr+".png")

def preprocess_dataframe():

    """
    Preprocesses the input data by normalizing the x1, x2, x3 values, computing the weighted values,
    and binning the data into a histogram.

    Parameters:
    -----------
    filename : str
        The name of the input CSV file, without the file extension.
    NUMBER_OF_NBINS : int
        The number of bins to use for the histogram calculation.

    Returns:
    --------
    df : pandas.DataFrame
        The preprocessed dataframe containing the normalized and binned data.
    binx : numpy.ndarray
        A linearly spaced array of bin edges in the range [0, 1].
    H : numpy.ndarray
        A histogram of the binned data, normalized by the sum of the histogram values.
    """

    df = pd.read_csv(filename+".csv", header=0, 
                     names=['ID','Total','x1','x2','x3','weight'])

    df['x1'] = df['x1']/100
    df['x2'] = df['x2']/100
    df['x3'] = df['x3']/100

    df['on_x1bin'] = -1
    df['on_x2bin'] = -1
    df['on_x3bin'] = -1

    df['x1_wb'] = df['x1']*df['weight']
    df['x2_wb'] = df['x2']*df['weight']
    df['x3_wb'] = df['x3']*df['weight']

    H, b = np.histogramdd((df['x1'].to_numpy(), df['x2'].to_numpy(), df['x3'].to_numpy()),
                          bins=(NUMBER_OF_NBINS, NUMBER_OF_NBINS, NUMBER_OF_NBINS), range=((0, 1), (0, 1), (0, 1)), weights=df['weight'].to_numpy())
    H = H / np.sum(H)

    binx = np.linspace(0, 1, NUMBER_OF_NBINS)

    for a,b_axis in enumerate(b):
        for i,t_bin in enumerate(b_axis):
            if i < 11:
                if t_bin == 0:
                    df.loc[(df['x'+str(a+1)] >= t_bin) & (df['x'+str(a+1)] <= b[a][i+1]), ['on_x'+str(a+1)+'bin']] = i
                else:
                    df.loc[(df['x'+str(a+1)] > t_bin) & (df['x'+str(a+1)] <= b[a][i+1]), ['on_x'+str(a+1)+'bin']] = i

    return df,binx,H

def generate_dicts(H, binx, R_dict):

    """
    Generates two dictionaries: one for the filtered triangular graph data and the other for the
    radii of the average METs at each point in the triangular graph.

    Parameters:
    -----------
    H : numpy.ndarray
        A histogram of the binned data, normalized by the sum of the histogram values.
    binx : numpy.ndarray
        A linearly spaced array of bin edges in the range [0, 1].
    R_dict : dict
        A dictionary of the radii for the average METs at each point in the triangular graph.

    Returns:
    --------
    interp_dict : dict
        A dictionary of the filtered triangular graph data. The data is filtered using a Gaussian
        filter if the 'args.filter' flag is set.
    bubble_dict : dict
        A dictionary of the radii of the average METs at each point in the triangular graph.
    """

    interp_dict = dict()
    bubble_dict = dict()

    kde = gaussian_filter(H, sigma=GAUSSIAN_FILTER_SIGMA)
    for i, x in enumerate(binx):
        for j, y in enumerate(binx):
            for k, z in enumerate(binx):

                if args.filter:
                    interp_dict[(i, j, k)] = kde[i, j, k]
                else:
                    interp_dict[(i, j, k)] = H[i, j, k]

                if i+j+k == 10:
                    bubble_dict[(i, j, k)] = R_dict[i, j, k]
                else:
                    bubble_dict[(i, j, k)] = 0

    return interp_dict, bubble_dict

def main():

    """
    Main function that executes the entire process for generating and exporting a triangular graph.

    The process includes:
    1. Preprocessing the input dataframe.
    2. Calculating the scatter radii for the triangular graph.
    3. Generating the dictionaries for the filtered triangular graph data and the bubble radii.
    4. Exporting the triangular graph with the processed data.
    """

    # Step 1: Preprocess the input dataframe
    df, binx, H = preprocess_dataframe()

    # Step 2: Calculate the scatter radii for the triangular graph
    R_dict = calculate_scatter_radius(df, binx)

    # Step 3: Generate the dictionaries for the filtered triangular graph data and the bubble radii
    interp_dict, bubble_dict = generate_dicts(H, binx, R_dict)

    # Step 4: Export the triangular graph with the processed data
    export_graph(df, interp_dict, bubble_dict)


if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument("--country", help="3 letter", type=str)
    parser.add_argument('--average', action='store_true', help="show average point")
    parser.add_argument('--filter', action='store_true', help="adapt gaussian_filter to ternary graph")
    parser.add_argument('--mets', action='store_true', help="show reference size of mets")
    parser.add_argument('--output', action='store_true', help="output csv file")
    args = parser.parse_args()


    filename = "data_"+str(args.country)
    optionstr = ""
    if args.average:
        optionstr = optionstr + "_average"
    if args.filter:
        optionstr = optionstr + "_filter"
    if args.mets:
        optionstr = optionstr + "_mets"

    main()
