import concurrent.futures
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import pyBigWig
from sklearn.cluster import KMeans


class TypeAnalysis(Enum):
    """
    Enum representing the type of analysis.

    Attributes
    ----------
    STRANDED : str
        Stranded analysis type.
    UNSTRANDED : str
        Unstranded analysis type.
    """

    STRANDED = "transcription"
    UNSTRANDED = "chip"


@dataclass
class MatrixAnalysis:
    """
    A dataclass representing matrix analysis options.

    Attributes
    ----------
    type_analysis : Optional[TypeAnalysis], optional
        The type of analysis to be performed. Defaults to None.
    k_means : Optional[int], optional
        The number of clusters for k-means clustering. Defaults to None.
    region_kmeans : Optional[Tuple[int, int]], optional
        The range of regions for k-means clustering. Defaults to None.
    """

    type_analysis: Optional[TypeAnalysis] = None
    k_means: Optional[int] = None
    region_kmeans: Optional[Tuple[int, int]] = None


@dataclass
class Quantification:
    """
    A class used to represent quantification data.

    Attributes
    ----------
    regions : str
        The path to the regions file.
    bigwigs : List[str]
        The list of paths to the bigwig files.
    name : str
        The name of the quantification.
    matrix : np.ndarray, optional
        The quantification matrix. Defaults to None.
    total_matrix_counts : Optional[float], optional
        The total matrix counts. Defaults to None.
    normalized_matrix : np.ndarray, optional
        The normalized matrix. Defaults to None.
    matrix_analysis : MatrixAnalysis, optional
        The matrix analysis options. Defaults to an instance of MatrixAnalysis.

    Methods
    -------
    make_matrix_and_get_total_counts(total_rows: int) -> Tuple[np.ndarray, float]:
        Makes the quantification matrix and returns it along with the total counts.
    normalize_matrix(norm_factor: int) -> None:
        Normalizes the matrix using the given normalization factor.
    """

    regions: str
    bigwigs: List[str]
    name: str
    matrix: np.ndarray = field(default=None)
    total_matrix_counts: Optional[float] = field(default=None)
    normalized_matrix: np.ndarray = field(default=None)
    matrix_analysis: MatrixAnalysis = MatrixAnalysis()

    @property
    def data_type(self):
        """
        Returns the analysis type.

        Returns
        -------
        str
            The name of the type of analysis.
        """
        return self.type_analysis.value

    def _get_counts_from_bigwig(self, args) -> pd.DataFrame:
        """
        Helper method to get counts from the bigwig files.

        Parameters
        ----------
        args : Tuple
            The arguments containing information about the bigwig files.

        Returns
        -------
        pd.DataFrame
            The dataframe of counts.
        float
            The sum of counts.
        """
        forward_bw, data_type, sliced_df, reverse_bw, antisense = args
        per_row_matrix_list = []
        # iterate through bed file of genomic regions
        for idx, each_row in sliced_df.iterrows():
            # get the coordinates of the gene
            chromosome, left_coordinate, right_coordinate, strand = (
                str(each_row[0]),
                int(each_row[1]),
                int(each_row[2]),
                str(each_row[5]),
            )
            if strand == "+":
                # if interested in antisense strand, use RV bigwig
                if antisense:
                    bw = pyBigWig.open(reverse_bw)
                else:
                    # intantiate FW bigwig for positive strand genes
                    bw = pyBigWig.open(forward_bw)
                # get values for each base: chromosome,start,end
                try:
                    per_row_matrix = np.array(
                        bw.values(chromosome, left_coordinate, right_coordinate)
                    )
                    per_row_matrix[np.isnan(per_row_matrix)] = 0

                except RuntimeError:
                    print(
                        f"Interval out of bound for {chromosome}:{left_coordinate}-{right_coordinate}"
                    )
                    per_row_matrix = np.zeros(right_coordinate - left_coordinate)
                    pass

            elif strand == "-":
                # if interested in antisense strand, use RV bigwig
                if antisense:
                    bw = pyBigWig.open(forward_bw)
                else:
                    # if dealing with chip data, data is only present in the forward bigwig
                    if data_type == TypeAnalysis.UNSTRANDED:
                        bw = pyBigWig.open(forward_bw)
                    else:
                        # intantiate RV bigwig for negative strand genes
                        bw = pyBigWig.open(reverse_bw)
                # get values for each base: chromosome,start,end
                try:
                    per_row_matrix = np.array(
                        bw.values(chromosome, left_coordinate, right_coordinate)
                    )
                    per_row_matrix[np.isnan(per_row_matrix)] = 0

                except RuntimeError:
                    print(
                        f"Interval out of bound for {chromosome}:{left_coordinate}-{right_coordinate}"
                    )
                    per_row_matrix = np.zeros(right_coordinate - left_coordinate)
                    pass

                # Reverse the order of elements in the array along horizontal axis for negative strand genes
                per_row_matrix = np.flip(per_row_matrix, 0)
            # Get the absolute value in case RV strand bw have 'negate' values activated.
            per_row_matrix_list.append(abs(np.array(per_row_matrix)))

        sum_counts = np.round(sum([sum(x) for x in per_row_matrix_list]), decimals=2)

        return pd.DataFrame(per_row_matrix_list).fillna(0), sum_counts

    def make_matrix_and_get_total_counts(
        self,
        total_rows: int,
        antisense: bool,
    ) -> Tuple[np.ndarray, float]:
        """
        Makes the quantification matrix and returns it along with the total counts.

        Parameters
        ----------
        total_rows : int
            The total number of rows.

        Returns
        -------
        Tuple[np.ndarray, float]
            The quantification matrix and the total counts.
        """
        genomic_regions = pd.read_csv(self.regions, sep="\t", header=None)
        fbwigs = self.bigwigs[0]
        rbwigs = self.bigwigs[1]
        how_analyze = self.matrix_analysis.type_analysis

        with concurrent.futures.ProcessPoolExecutor() as executor:
            print("Making gene matrix...")
            results = list(
                executor.map(
                    self._get_counts_from_bigwig,
                    [
                        (
                            fbwigs,
                            how_analyze,
                            genomic_regions[rows : rows + 100],
                            rbwigs,
                            antisense,
                        )
                        for rows in range(0, total_rows + 100, 100)
                    ],
                )
            )

        # cocatenate row matrices and sum reads
        matrices_per_row = (pd.concat([x[0] for x in results])).to_numpy()
        sums_per_row = sum([x[1] for x in results])
        print("Done making matrix!")

        return matrices_per_row, sums_per_row

    def normalize_matrix(self, norm_factor: int) -> None:
        """
        Normalizes the matrix using the given normalization factor.

        Parameters
        ----------
        norm_factor : int
            The normalization factor.

        Returns
        -------
        None
        """
        print(
            f"The total number of reads in {self.name} is {self.total_matrix_counts:,.2f}. ",
            end="",
        )

        self.normalized_matrix = np.round(self.matrix * norm_factor, decimals=2)


def kmeans_clustering(
    bed_path,
    output_path,
    sliced_fc_matrix,
    num_matrix,
    deno_matrix,
    fc_matrix,
    k,
):
    """
    Performs k-means clustering on the matrices.

    Parameters
    ----------
    sliced_fc_matrix : np.ndarray
        The sliced fold-change matrix.
    num_matrix : np.ndarray
        The numerator matrix.
    deno_matrix : np.ndarray
        The denominator matrix.
    fc_matrix : np.ndarray
        The fold-change matrix.
    k : int
        The number of clusters.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        The reordered matrices after clustering.
    """
    print("\nClustering...")
    bed_file = pd.read_csv(bed_path, sep="\t", header=None)
    kmeans = KMeans(n_clusters=k).fit(sliced_fc_matrix)
    # order the centroids in an attempt to
    # get the same cluster order
    cluster_labels = kmeans.labels_

    # add the cluster labels to all matrices based on the control matrix
    clustered_matrix_labels_num = np.concatenate(
        (num_matrix, np.array([cluster_labels]).T), axis=1
    )
    clustered_matrix_labels_deno = np.concatenate(
        (deno_matrix, np.array([cluster_labels]).T), axis=1
    )
    clustered_matrix_labels_fc = np.concatenate(
        (fc_matrix, np.array([cluster_labels]).T), axis=1
    )
    clustered_bed_file = np.concatenate(
        (bed_file.to_numpy(), np.array([cluster_labels]).T), axis=1
    )

    # re-order the matrix based on the cluster labels
    re_ordered_matrix_num = clustered_matrix_labels_num[
        clustered_matrix_labels_num[:, -1].argsort()
    ]
    re_ordered_matrix_deno = clustered_matrix_labels_deno[
        clustered_matrix_labels_deno[:, -1].argsort()
    ]
    re_ordered_matrix_fc = clustered_matrix_labels_fc[
        clustered_matrix_labels_fc[:, -1].argsort()
    ]

    re_ordered_bed_file = clustered_bed_file[clustered_bed_file[:, -1].argsort()]

    # remove last column (cluster labels)
    re_ordered_matrix_num = re_ordered_matrix_num[:, :-1]
    re_ordered_matrix_deno = re_ordered_matrix_deno[:, :-1]
    re_ordered_matrix_fc = re_ordered_matrix_fc[:, :-1]

    # output as pandas dataframe
    re_ordered_bed_file = pd.DataFrame(re_ordered_bed_file)
    re_ordered_bed_file.to_csv(
        f"{output_path}Clustered_genes.bed", sep="\t", index=False, header=False
    )

    return re_ordered_matrix_num, re_ordered_matrix_deno, re_ordered_matrix_fc
