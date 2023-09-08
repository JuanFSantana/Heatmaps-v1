import argparse
import concurrent.futures
import sys
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import pyBigWig

from utils.check_bed_utils import BedFile
from utils.heatmap_utils import Heatmap, HeatmapAnalysisType, HeatmapColorType
from utils.quantify_utils import (
    TypeAnalysis,
    MatrixAnalysis,
    Quantification,
    kmeans_clustering,
)


def imaging(heatmap_obj: Heatmap):
    """
    Creates a heatmap image using the input Heatmap object's image method.

    Parameters:
    heatmap_obj (Heatmap): An instance of a Heatmap object.

    Returns:
    Object: The image object created by the heatmap's image method.
    """
    print(f"Creating heatmaps for {heatmap_obj.identifier}...")
    return heatmap_obj.image()


def making_images(heatmap_obj: Heatmap) -> None:
    """
    Concurrently creates images for each Heatmap object.

    Parameters:
    heatmap_obj (Heatmap): An instance of a Heatmap object.

    Returns:
    None
    """
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(imaging, heatmap_obj)


def modifiy_base_per_pixel(otro_tres_tuple):
    """
    Modifies a matrix based on the width and height provided. It supports repeating or averaging the matrix
    both vertically and horizontally.

    Parameters:
    otro_tres_tuple (Tuple): A tuple containing a DataFrame, height, and width.

    Returns:
    np.ndarray: The modified matrix as a numpy array.
    """
    table, height, width = otro_tres_tuple

    # Aspect of height and width
    if height >= 1 and width >= 1:
        verically_repeated = table.reindex(table.index.repeat(height))
        horizontally_repeated = verically_repeated.reindex(
            verically_repeated.columns.repeat(width), axis="columns"
        )
        final_matrix = horizontally_repeated

    elif height >= 1 and width < 1:
        # average first
        # pixels per base
        width_rolling_avg = int(1 / width)
        # rolling average and select rows containing the average window HORIZONTALLY
        df_matrix_width_avg = (
            table.rolling(width_rolling_avg, axis=1).mean().dropna(axis=1, how="any")
        )
        avg_matrix = df_matrix_width_avg[
            df_matrix_width_avg.columns[::width_rolling_avg]
        ]
        # repeat array vertically
        verically_repeated = avg_matrix.reindex(avg_matrix.index.repeat(height))
        final_matrix = verically_repeated

    elif height < 1 and width < 1:
        height_rolling_avg = int(1 / height)
        width_rolling_avg = int(1 / width)
        # rolling average and select rows containing the average window VERTICALLY
        height_avg_matrix = (
            table.rolling(height_rolling_avg, axis=0)
            .mean()
            .dropna(axis=0, how="any")
            .iloc[::height_rolling_avg]
        )
        # rolling average and select rows containing the average window HORIZONTALLY
        height_width_avg_matrix = (
            height_avg_matrix.rolling(width_rolling_avg, axis=1)
            .mean()
            .dropna(axis=1, how="any")
        )
        avg_matrix = height_width_avg_matrix[
            height_width_avg_matrix.columns[::width_rolling_avg]
        ]
        final_matrix = avg_matrix

    elif height < 1 and width >= 1:
        # average first
        # pixels per base
        height_rolling_avg = int(1 / height)
        # rolling average and select rows containing the average window VERTICALLY
        height_avg_matrix = (
            table.rolling(height_rolling_avg, axis=0)
            .mean()
            .dropna(axis=0, how="any")
            .iloc[::height_rolling_avg]
        )
        # repeat array horizontally
        horizontally_repeated = height_avg_matrix.reindex(
            height_avg_matrix.columns.repeat(width), axis="columns"
        )
        final_matrix = horizontally_repeated
    print("Averaging/repeating matrix")
    return final_matrix.to_numpy()


def matrix_division(numerator, denominator):
    """
    Calculates the log2 fold change of a division operation between two arrays.

    Parameters:
    numerator (np.ndarray): The numerator array.
    denominator (np.ndarray): The denominator array.

    Returns:
    np.ndarray: Array representing log2 fold change.
    """
    fold_change_array = np.divide(
        numerator,
        denominator,
        out=np.array(numerator),
        where=denominator != 0,
        dtype=float,
    )
    division_log2_array = np.log2(
        fold_change_array,
        out=np.array(fold_change_array),
        where=fold_change_array != 0,
        dtype=float,
    )
    print("\nDone calculating log2 fold change")

    return division_log2_array


def conccurent_matrix_maker(args) -> None:
    """
    Concurrently creates matrices using the input Quantification object and the specified number of rows.

    Parameters:
    args (Tuple): A tuple containing a Quantification object and the number of rows.

    Returns:
    None
    """
    quant_obj, num_rows, antisense = args
    return quant_obj.make_matrix_and_get_total_counts(num_rows, antisense)


def make_matrix(
    nume_obj: Quantification, deno_obj: Quantification, num_rows: int, anti_sense: bool
) -> None:
    """
    Makes two matrices concurrently using the specified Quantification objects and the number of rows.
    Then, it assigns the resulting matrices and total counts to the respective Quantification objects.

    Parameters:
    nume_obj (Quantification): The Quantification object for the numerator.
    deno_obj (Quantification): The Quantification object for the denominator.
    num_rows (int): The number of rows for the matrices.
    anti_sense (bool): If True, the anti-sense strand will be used instead of the sense strand.

    Returns:
    None
    """
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = list(
            executor.map(
                conccurent_matrix_maker,
                [(obj, num_rows, anti_sense) for obj in [nume_obj, deno_obj]],
            )
        )
    # results is a list of tuples, where each tuple contains a matrix and total counts
    nume_obj.matrix, nume_obj.total_matrix_counts = results[0]
    deno_obj.matrix, deno_obj.total_matrix_counts = results[1]


def parse_args():
    """
    Get arguments
    """

    parser = argparse.ArgumentParser(
        prog="create_heatmap.py",
        description="Generates heatmaps of counts per base over a chosen genomic interval for control, experimental and log2FC",
    )
    parser.add_argument(
        "regions",
        type=str,
        help="Bed file of genomic regions of chosen length. Gene/feature name or other identifier must be in the 4th column",
    )
    parser.add_argument(
        "-numerator",
        dest="numerator",
        metavar="\b",
        required=True,
        type=str,
        nargs="*",
        help="Forward and reverse (in that order) bigwigs for condition to be used as numerator in log2FC",
    )
    parser.add_argument(
        "-denominator",
        dest="denominator",
        metavar="\b",
        required=True,
        type=str,
        nargs="*",
        help="Forward and reverse (in that order) bigwigs for condition to be used as denominator in log2FC",
    )
    parser.add_argument(
        "-mb",
        dest="max_black_value",
        metavar="\b",
        default="default",
        nargs="*",
        help="Sets the chosen value as black, default is largest number in the matrix",
    )
    parser.add_argument(
        "-mc",
        dest="max_color_value",
        metavar="\b",
        default="default",
        nargs="*",
        help="Sets the chosen value as red/blue, default is absolute largest number in the matrix",
    )
    parser.add_argument(
        "-chip",
        dest="chip_data",
        action="store_true",
        default=False,
        help="If bigwigs are ChIP-seq data, argument must be invoked",
    )
    parser.add_argument(
        "-k",
        dest="kmeans",
        metavar="\b",
        default=None,
        nargs=1,
        type=int,
        help="Number of K-means clusters to be used for clustering the data. A .csv file with the cluster assignments will be created in the output directory.",
    )

    parser.add_argument(
        "-r",
        dest="region_kmean",
        metavar="\b",
        nargs=2,
        type=int,
        help="Location in the regions file to be used for k-means clustering. The first argument is the start position and the second is the end position.\
                            For example, in regions of 1000 bp, if you want to cluster the middle 500 bp, the arguments would be -r 250 750. ",
    )

    parser.add_argument(
        "-y",
        dest="y_axis",
        metavar="\b",
        type=float,
        default=1.0,
        help="Horizontal lines/bp for each row displayed, for example -y 1 is one horizontal line/region or gene; -y 0.1 is one averaged horizontal line/10 region or gene",
    )
    parser.add_argument(
        "-x",
        dest="x_axis",
        metavar="\b",
        type=float,
        default=1.0,
        help="Vertical lines/bp for each row displayed, for example -y 1 is one horizontal line/bp; -y 0.1 is one averaged horizontal line/10 bp",
    )
    parser.add_argument(
        "-g",
        dest="gamma",
        metavar="\b",
        type=float,
        default=1.0,
        help="Gamma correction",
    )
    parser.add_argument(
        "-n",
        dest="normalize",
        action="store_true",
        default=False,
        help="If argument is invoked, the average total number of reads for all regions will be calculated between the numerator and denominator datasets and the reads per base will be normalized to this value",
    )

    parser.add_argument(
        "-antiSense",
        dest="antiSense",
        action="store_true",
        default=False,
        help="If argument is invoked, the counts from the anti-sense strand will be plotted instead of the sense strand. This is only applicable for transcriptional data.",
    )

    parser.add_argument(
        "-o",
        dest="output_dir",
        metavar="\b",
        type=str,
        required=True,
        nargs=1,
        help="Path to output",
    )

    parser.add_argument(
        "-id",
        dest="names",
        metavar="\b",
        type=str,
        required=True,
        nargs=2,
        help="Image output names for numerator and denominator (in that order). The log2FC heatmaps will be named by combining the numerator and denominator file names",
    )

    args = parser.parse_args()

    file_regions = args.regions
    max_val = args.max_black_value
    max_color = args.max_color_value
    height = float(args.y_axis)
    width = float(args.x_axis)
    gamma = float(args.gamma)
    kmeans = int(args.kmeans[0]) if args.kmeans else None
    region_kmean = args.region_kmean
    to_normalize = args.normalize
    chip_data = args.chip_data
    output_directory = args.output_dir[0]
    numerator_id, denominator_id = args.names
    anti_sense = args.antiSense

    if chip_data:
        if anti_sense:
            sys.exit(
                "The -antiSense argument is only applicable for transcriptional data."
            )
        if len(args.numerator) != 1 or len(args.denominator) != 1:
            sys.exit(
                "For ChIP-seq data, only one bigwig file is needed for each condition in arguments -numerator and -denominator. If transcriptional data is used, two bigwig files are needed for each condition."
            )
        else:
            bw_control_fw, bw_control_rv = args.numerator[0], None
            bw_experimental_fw, bw_experimental_rv = args.denominator[0], None
    else:
        if len(args.numerator) != 2 or len(args.denominator) != 2:
            sys.exit(
                "For transcriptional data, two bigwig files are needed for each condition in arguments -numerator and -denominator. If ChIP-seq data is used, invoke the -chip argument, and only one bigwig file is needed for each condition."
            )
        else:
            bw_control_fw, bw_control_rv = args.numerator
            bw_experimental_fw, bw_experimental_rv = args.denominator

    start_kmean, end_kmean = None, None
    if kmeans:
        if not region_kmean:
            sys.exit(
                "If k-means clustering is to be used, the -r argument must be invoked."
            )
        else:
            start_kmean, end_kmean = region_kmean
            start_kmean, end_kmean = int(start_kmean), int(end_kmean)

    if region_kmean:
        if not kmeans:
            sys.exit(
                "If the -r argument is invoked, the -k argument must also be invoked."
            )

    args = [
        file_regions,
        bw_control_fw,
        bw_control_rv,
        bw_experimental_fw,
        bw_experimental_rv,
        max_val,
        height,
        width,
        gamma,
        output_directory,
        to_normalize,
        chip_data,
        max_color,
        kmeans,
        start_kmean,
        end_kmean,
        numerator_id,
        denominator_id,
        anti_sense,
    ]

    return args


def main(args):
    (
        file_regions,
        bw_numerator_fw,
        bw_numerator_rv,
        bw_denominator_fw,
        bw_denominator_rv,
        max_val,
        height,
        width,
        gamma,
        output_directory,
        to_normalize,
        chip_data,
        max_color,
        kmeans,
        start_kmean,
        end_kmean,
        numerator_id,
        denominator_id,
        anti_sense,
    ) = args

    # check bed file
    num_rows, num_cols = BedFile(file_regions).check_regions()

    # Creating Quantification objects
    # create bigwig quantification objects
    nume_bw_obj = Quantification(
        regions=file_regions,
        bigwigs=[bw_numerator_fw, bw_numerator_rv],
        name=numerator_id,
    )

    deno_bw_obj = Quantification(
        regions=file_regions,
        bigwigs=[bw_denominator_fw, bw_denominator_rv],
        name=denominator_id,
    )

    # Creating a MatrixAnalysis object
    shared_matrix_analysis = MatrixAnalysis(
        k_means=kmeans,
        region_kmeans=(start_kmean, end_kmean),
        type_analysis=TypeAnalysis.UNSTRANDED if chip_data else TypeAnalysis.STRANDED,
    )
    # Update the shared matrix_analysis for all Quantification instances
    nume_bw_obj.matrix_analysis = shared_matrix_analysis
    deno_bw_obj.matrix_analysis = shared_matrix_analysis

    # make matrix
    make_matrix(nume_bw_obj, deno_bw_obj, num_rows, anti_sense)
    # normalize matrices
    if to_normalize:
        # calculate normalization factor
        avg_counts = (
            nume_bw_obj.total_matrix_counts + deno_bw_obj.total_matrix_counts
        ) / 2
        # normalize matrices
        for obj in [nume_bw_obj, deno_bw_obj]:
            normalization_factor = avg_counts / obj.total_matrix_counts
            obj.normalize_matrix(normalization_factor)
            print(f"The normalization factor is {normalization_factor:,.2f}")

        # calculate log2FC
        numerator_matrix = nume_bw_obj.normalized_matrix
        denominator_matrix = deno_bw_obj.normalized_matrix

    else:
        # calculate log2FC
        numerator_matrix = nume_bw_obj.matrix
        denominator_matrix = deno_bw_obj.matrix

    # calculate log2FC
    log2FC = matrix_division(numerator_matrix, denominator_matrix)
    # k-means clustering
    if shared_matrix_analysis.k_means:
        sliced_fc_matrix = log2FC[
            :,
            shared_matrix_analysis.region_kmeans[
                0
            ] : shared_matrix_analysis.region_kmeans[1],
        ]
        print(
            f"\nKmeans clustering analysis is being applied to positions {shared_matrix_analysis.region_kmeans[0]} to {shared_matrix_analysis.region_kmeans[1]} of the log2FC matrix"
        )
        nume_clustered, deno_clustered, fc_clustered = kmeans_clustering(
            file_regions,
            output_directory,
            sliced_fc_matrix,
            numerator_matrix,
            denominator_matrix,
            log2FC,
            shared_matrix_analysis.k_means,
        )

    # calculate averages/ repeat pixels
    matrices = [
        pd.DataFrame(
            nume_clustered if shared_matrix_analysis.k_means else numerator_matrix
        ),
        pd.DataFrame(
            deno_clustered if shared_matrix_analysis.k_means else denominator_matrix
        ),
        pd.DataFrame(fc_clustered if shared_matrix_analysis.k_means else log2FC),
    ]
    with concurrent.futures.ProcessPoolExecutor() as executor:
        mod_matrices = list(
            executor.map(
                modifiy_base_per_pixel,
                [(each_matrix, height, width) for each_matrix in matrices],
            )
        )

    # create heatmap objects
    heatmap = []
    names = [numerator_id, denominator_id, f"{numerator_id}_{denominator_id}"]
    for idx, array in enumerate(mod_matrices):
        if idx == 0 or idx == 1:
            if chip_data:
                color_type = HeatmapColorType.CHIP
            else:
                color_type = HeatmapColorType.BLACKNWHITE
            name = names[idx]
        elif idx == 2:
            color_type = HeatmapColorType.FOLDCHANGE
            name = names[idx]
        heatmap.append(
            Heatmap(
                array=array,
                heatmap_type=color_type,
                heatmap_analysis_type=HeatmapAnalysisType.CLASSIC,
                max_value=max_val,
                max_color_value=max_color,
                identifier=name,
                y_axis=height,
                x_axis=width,
                yaxis_min=0,
                yaxis_max=num_rows,
                gamma=gamma,
                output_directory=output_directory,
            )
        )

    # if the max values is set to 'default', use the the largest max value between the numerator and denominator and set that list as the max value for both objects
    if isinstance(heatmap[0].max_value, list):
        num_max_val = heatmap[0].max_value[0]
        deno_max_val = heatmap[1].max_value[0]

        if num_max_val > deno_max_val:
            heatmap[1].max_value = heatmap[0].max_value
        else:
            heatmap[0].max_value = heatmap[1].max_value

    # create heatmaps
    making_images(heatmap)


if __name__ == "__main__":
    args = parse_args()
    main(args)
