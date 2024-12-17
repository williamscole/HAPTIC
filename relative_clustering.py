# from concurrent.futures import process
import concurrent.futures
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import networkx as nx
import itertools as it
import sys
from collections import Counter
import time
from helper_functions import annotate_haplotypes, segment_coverage, get_roh_from_relative, column_swapper, split_regions
import math
import pickle as pkl
import warnings 

warnings.simplefilter(action='ignore', category=FutureWarning)

#last update: 9fdccba
# before 7JUN roh update: 4112533
# 

def parse_args():
    """Parses the command line arguments and returns the ArgumentParser.parse_args() object."""
    parser = argparse.ArgumentParser()
    ### inputs

    # only for use if there are batches
    # otherwise has the following format: [batch keyword],[cur batch (batch where the current focal indv are from)],[range of batch numvbers]
    # e.g., batch,1,1-26 means that our focal individuals are in batch1 and we're loading batches 1 through 26
    # e.g., batch,1,1-10,12-20,26 would mean that we're loading batches 1 thru 10, 12 thru 20 and 26
    parser.add_argument("-batches", type=str, default = ",,")

    parser.add_argument("-focal_file", help = "First column is IDs of the individuals in the run; 2nd and 3rd cols are optional, but are for listing the parents")

    # path file containing the IBD segments
    parser.add_argument("-relative_ibd_file")

    # path file containing the IBD segments shared with the parent
    parser.add_argument("-parent_ibd_file", default = "")

    # if all chrom are in one file (doesn't look for chr1, chr2, etc)
    parser.add_argument("-single_ibd_file", help = "For use if all chromosomes present in one IBD file", action = "store_true")

    ### various params
    parser.add_argument("-min_seg_length", help="minimum cMs segment length for inclusion", type=float, default=5.0)
    parser.add_argument("-min_k", help="Minimum kinship in cM with focal for inclusion", type=float, default=10)
    parser.add_argument("-max_k", help="Maximum kinship in cM with focal for inclusion", type=float, default=2000)


    ### ignores tiles on the X chrom
    parser.add_argument("-noX", action = "store_true")

    ### what to name the x
    parser.add_argument("-x_name", type=str, default="23")

    ### output results prefix
    parser.add_argument("-out", type=str, default="clustering")

    ### don't remove ROH; don't remove high IBD density regions
    parser.add_argument("-keepROH", action = "store_true")

    # standard devs 
    parser.add_argument("-stdevs", help="For IBD segments, computes mean IBD coverage and removes IBD in high-density regions, defined as being [stdevs] sds away from mean.", default=3, type=int)

    ### compute ROH from the relatives
    parser.add_argument("-relative_roh", action = "store_true")

    ### sex file
    parser.add_argument("-sex_file", type = str, default = None)

    ## turn developer mode on
    parser.add_argument("-dev_mode", action = "store_true")

    ### debug mode skips the multiprocessing
    parser.add_argument("-multiprocess", action = "store_true")

    ### only write out the IBD segments (useful for parallelizing the clustering runs)
    parser.add_argument("-ibd_only", action = "store_true")

    ### can supply the pickled IBD file
    parser.add_argument("-ibd_pkl", default="")

    ### maximum kinship for the descendant-free runs
    parser.add_argument("-descendant_k", help="For the initial run to find descendants, the max-k of relatives included. Set to 0 to skip the descendant-finder.", default=1100.0, type=float)

    ### use if only want the haplotype annotation
    parser.add_argument("-hap_annotation", action="store_true")

    ### choose the num of threads, which is to say the size of the chunks to break focals in
    parser.add_argument("-threads", type=int, default=5)

    ### should be a speed up for 23andMe; read the help
    parser.add_argument("-coverage_gain", help="When iterating through different thresholds, if the smallest eigenvalue is >0 and the gain in coverage is less than this, then stop clustering, i.e., if (cur_cov-prev_cov)/cur_cov < coverage_gain.",
                        type=float, default=-np.inf)
    
    ### if you want to continue the run
    parser.add_argument("-continue_run", help="If already run, will continue writing to {out}_results.pkl", action="store_true")

    ### customize the t2 thesholds
    parser.add_argument("-t2", nargs="*", type=int, default=[10, 25, 50, 75, 100])

    args = parser.parse_args()

    return args


'''IBDHandler loads IBD segments from phasedibd feather files.
It removes small segments, close or distant relatives, removes ROH
and high density IBD regions'''
class IBDHandler:

    def __init__(self, args):

        # holds all the arguments
        self.args = args

        # empty roh_df
        self.roh = {F: pd.DataFrame(columns=["id1", "chromosome", "start_cm", "end_cm"])\
                                    for F in args.focal_ids}

        # empty high coverage regions
        self.high_cov = pd.DataFrame(columns=["start_cm", "end_cm", "chromosome"])

        # store focal ibd
        self.ibd = {F: pd.DataFrame(columns=["id1", "id2", "id1_haplotype", "id2_haplotype",
                                "start_cm", "end_cm", "l", "chromosome"]) for F in args.focal_ids}

    ### loads pandas feather files, ensures same dtype, sorts id1, id2 columns appropriately
    ### doesn't do much with data: excludes non-focal samples, small segments
    ### also stores the data in the IBDHandler class
    def load_ibd_file(self, filename, multiprocess=False):

        # default is to call the chromosome 23, but can replace it with X
        filename = filename.replace("chr23", f"chr{self.args.x_name}")

        print(filename)

        # read in ibd and reset the index
        df = pd.read_feather(filename)

        # get the segment lengths and remove small segments
        df["l"] = df["end_cm"] - df["start_cm"]
        df = df[df.l > self.args.min_seg_length].reset_index(drop=True)

        ### ensure chrom is an int
        if type(df["chromosome"][0]) == str:

            df["chromosome"] = df["chromosome"].apply(lambda x: x[3:] if "chr" in x else x)

        df["chromosome"] = df["chromosome"].apply(lambda x: int(x) if x != "X" else 23)

        ### ensure ids are str
        df["id1"] = df["id1"].apply(str)
        df["id2"] = df["id2"].apply(str)

        # make sure that (1) only focal IDs' IBD in there and (2) id1 is the focal id
        tmp1 = df[df.id1.isin(self.args.focal_ids)].copy()
        tmp2 = df[df.id2.isin(self.args.focal_ids)].copy()

        # switch the columns for tmp2
        tmp2[["id1", "id2", "id1_haplotype", "id2_haplotype"]] = tmp2[["id2", "id1", "id2_haplotype", "id1_haplotype"]]

        # concat the two, remove duplicate rows (e.g., if id1 and id2 are both focal ids)
        df = pd.concat([tmp1, tmp2]).drop_duplicates(keep="first", ignore_index=True)

        # add each focal ibd
        if multiprocess:
            return df

        for focal, focal_df in df.groupby("id1"):
            # get the current ibd dataframe
            tmp = pd.concat([self.ibd[focal], focal_df]).reset_index(drop=True)
            # add the new to the current and set chromosome to int
            tmp["chromosome"] = tmp["chromosome"].apply(int)
            self.ibd[focal] = tmp

    def load_ibd_from_list(self, ibd_file_list, multiprocess=False):

        with concurrent.futures.ProcessPoolExecutor() as executor:
            # multiprocess loading the files
            if multiprocess:
                results = [executor.submit(self.load_ibd_file, ibd_file, True) for ibd_file in ibd_file_list]
                
                for r in concurrent.futures.as_completed(results):
                    # get the segments and add
                    for focal, focal_df in r.result().groupby("id1"):
                        # get previous focal df and concat
                        tmp = pd.concat([self.ibd[focal], focal_df]).reset_index(drop=True)
                        # make sure chromosome type is int
                        tmp["chromosome"] = tmp["chromosome"].apply(int)
                        # reassign the segments
                        self.ibd[focal] = tmp

            # no multiprocessing; run one file at a time
            else:
                for ibd_file in ibd_file_list:
                    self.load_ibd_file(ibd_file, multiprocess=False)


    '''
    Given IBD segments and regions to remove for a specific chromosome,
    removes those regions from the IBD segments and discards them if they
    are now too small
    '''
    def remove_regions_from_chromosome(self, segments_df, regions_df):
            
        # keeps index of segments that overlap with the region
        r_segment_idx = set()

        # iterate through the roh
        for _, region in regions_df.iterrows():
            # get start, stop pos of the region
            r1, r2 = region["start_cm"], region["end_cm"]

            # subset segments that overlap
            overlapped = segments_df[segments_df.apply(lambda x: min(r2, x.end_cm) - max(r1, x.start_cm) > 0, axis=1)].index
            r_segment_idx |= set(overlapped)

        # get the segments with no ROH; this will be the output
        no_region_chrom = segments_df[~segments_df.index.isin(r_segment_idx)]

        # dict for the split regions function
        region_d = {(row["start_cm"], row["end_cm"]): ["region"] for _, row in regions_df.iterrows()}

        # iterate through the segments that overlap with regions
        for _, segment in segments_df.loc[list(r_segment_idx)].iterrows():

            # split the segments across regions
            tmp = split_regions(region_d, [segment["start_cm"], segment["end_cm"], "ibd"])

            # get regions that are (1) only IBD and (2) long enough
            regions = [(start, stop) for (start, stop), rtype in tmp.items()\
                        if rtype == ["ibd"] and stop - start >= self.args.min_seg_length]

            # no segment remains
            if len(regions) == 0:
                continue

            # create new df and add the new start/stop
            tmp_df = pd.DataFrame([segment for _ in regions])
            tmp_df[["start_cm", "end_cm"]] = regions

            # add to the no region
            no_region_chrom = pd.concat([no_region_chrom, tmp_df])

        return no_region_chrom.reset_index(drop=True)


    ### given a df of IBD segments, does several things
    ### must be all chromosomes but could be just a single batch
    ### 1) removes small segments
    ### 2) computes kinship coefficients and removes individuals <k>
    ### 3) renames ibd segments
    def process_relative_ibd(self, focal, focal_segments, min_k, max_k):

        # takes as input the id name, the list of segments, and the existing index to id dict
        def add_name(id2, k, segments, index_to_id, id_n_segments):

            # n keeps track of the segment number
            n = id_n_segments.get(id2, 1)

            '''segment_coverage merges IBD segments that are separated by at most 1 cM gap
            It doesnt actually merge them in the IBD data frame; the merge is reflected in the ID of the segment
            Say you have the two following segments: [0, 10] and [10, 20] shared with relative A
            We would keep the two segments but they would both be named (A, 1, k, 20). The merging here is just telling us we 
            are very confident that these IBD segments were inherited from the same parent'''
            for start, stop, index in segment_coverage(segments, 1):

                # dont care about precision, just taking the ind
                l = int(stop - start)

                # index gives a list of indices associated with the merged segment
                for i in index:
                    # a segment is given a name which is a tuple of the id name, the segment number, kinship in cM, and the length of the merged seg
                    index_to_id[i] = (id2, n, k, l)
                # add another segment
                n += 1
            # update the dict
            id_n_segments[id2] = n

            return index_to_id, id_n_segments
        
        
        # no segments
        if focal_segments.shape[0] == 0:
            return focal, focal_segments[["id2","id1_haplotype","chromosome","start_cm","end_cm"]]

        # create dict where id maps to the kinship coeff in cM
        k_dict = {id2: int(id2_df["l"].sum()) for id2, id2_df in focal_segments.groupby("id2")}

        ### only keep individulas who are > min_k but less than max_k
        temp = focal_segments[focal_segments.apply(lambda  x: min_k < k_dict[x.id2] < max_k, axis=1)].reset_index(drop=True)

        if temp.shape[0] == 0:
            return focal, temp

        ### drop ROH and high_coverage regions

        # remove the roh and high coverage regions
        remove_df = pd.concat([self.roh[focal], self.high_cov])

        # remove the regions from each chromosome and concat the dfs
        temp = pd.concat([self.remove_regions_from_chromosome(chrom_df, remove_df[remove_df.chromosome==chrom])\
                          for chrom, chrom_df in temp.groupby("chromosome")])

        # sort by position
        temp = temp.sort_values(by="end_cm").sort_values(by="start_cm").reset_index(drop=True)

        ### rename each segment
        temp["index"] = temp.index

        # maps the index of the segment to the new ID name
        index_to_id = {}
        # for each id, keeps track of the number of segments
        id_n_segments = {}

        # iterate through segments
        for (id2, _), segment_df in temp.groupby(["id2","chromosome"]):

            # update the index to id dict one id2 at a time
            index_to_id, id_n_segments = add_name(id2, k_dict[id2], segment_df[["start_cm", "end_cm", "index"]].values.tolist(), index_to_id, id_n_segments)

        # convert the old name to the new name
        temp["id2"] = temp["index"].apply(lambda x: index_to_id[x])

        return focal, temp[["id2","id1_haplotype","chromosome","start_cm","end_cm"]]

    '''Computes regions of the genome with high IBD coverage
    and removes them. Removes focal-specific ROH before computing.'''
    def set_high_coverage_regions(self):

        '''given a focal individual, returns a matrix with the counts
        of IBD segments that span each cM integer'''
        def focal_cm_ints(focal):

            # get the segments and roh
            focal_segments = self.ibd[focal]
            focal_roh = self.roh[focal]

            # holds the counts for each site at each chromosome
            chrom_counts = np.zeros((290, 23))

            for chrom, chrom_segments in focal_segments.groupby("chromosome"):

                chrom_roh = focal_roh[focal_roh.chromosome==chrom]
                
                # get the integers of the the cMs covered by IBD
                ints = it.chain(*chrom_segments.apply(lambda x: np.arange(math.ceil(x.start_cm), math.floor(x.end_cm)), axis=1))

                # get the integers of the cMs with roh
                roh_ints = it.chain(*chrom_roh.apply(lambda x: np.arange(math.ceil(x.start_cm), math.floor(x.end_cm)), axis=1))\
                            if chrom_roh.shape[0] > 0 else []

                for i, i_count in Counter(ints).items():
                    # integer in an roh; don't count it
                    if i in roh_ints:
                        continue

                    # fill in the value in the matrix
                    chrom_counts[i, chrom-1] = i_count
            
            return chrom_counts

        # chrom_counts
        chrom_counts = sum([focal_cm_ints(focal) for focal in self.ibd])

        # get means, std, upper limit for coverage
        mean_count = np.mean(chrom_counts, where=chrom_counts > 0, axis=0)
        std_count = np.std(chrom_counts, where=chrom_counts > 0, axis=0)
        upper = mean_count + std_count*self.args.stdevs

        # np.where returns the sites and the chromosomes where there is high IBD coverage
        sites, chroms = np.where(chrom_counts > upper)

        # create df of the sites
        high_cov = pd.DataFrame(np.array([sites, sites+1, chroms+1]).T, columns=["start_cm", "end_cm", "chromosome"])

        # create a df of the regions
        high_cov_regions = pd.DataFrame(it.chain(*[segment_coverage(df.values.tolist(), 1e-10)\
                                                    for _,df in high_cov.groupby("chromosome")]),
                                                    columns=["start_cm", "end_cm", "chromosome"])
        # rename the chromosome column
        high_cov_regions["chromosome"] = high_cov_regions["chromosome"].apply(lambda x: list(x)[0])
        
        # assign to the class
        self.high_cov = high_cov_regions

    def process_all_ibd(self, multiprocess=False):

        # add the roh
        if not self.args.keepROH:
            for focal, focal_segments in self.ibd.items():
                # for each focal, compute and add the roh
                self.roh[focal] = get_roh_from_relative(focal_segments, overlap_cm=3)

        # get the high coverage regions
        self.set_high_coverage_regions()

        with concurrent.futures.ProcessPoolExecutor() as executor:
            # multiprocess
            if multiprocess:

                # multiprocess but don't overwhelm
                flist = list(self.ibd)
                for i in range(0, len(flist), 23):
                    print(f"Processing IBD segments ({i} of {len(flist)}")
                    results = [executor.submit(self.process_relative_ibd, focal, self.ibd[focal], self.args.min_k, self.args.max_k)\
                            for focal in flist[i:i+23]]
                    
                    for r in concurrent.futures.as_completed(results):
                        # get the segments and add
                        focal, focal_df = r.result()

                        if focal_df.shape[0] > 0:
                            self.ibd[focal] = focal_df

            # no multiprocessing; run one focal at a time
            else:

                for focal, focal_segments in self.ibd.items():
                    # get the segments
                    _, focal_df = self.process_relative_ibd(focal, focal_segments, self.args.min_k, self.args.max_k)
                    self.ibd[focal] = focal_df





'''For a focal individual, this class stores and computes tiles. It has a major upgrade over
the previous version, which could not be updated; e.g., in this version you can update the tiles
whenever and create a new tile_df. Note that you can add more segments to it, but you cannot take segments
away.'''
class Tiles:

    def __init__(self, cm_gap=3):

        # a list for each chrom
        self.chrom = {chrom: [] for chrom in range(1, 24)}

        # maps the index of an ibd segment to the relative ID
        self.index_to_id = {}

        self.gap = cm_gap

    def add_chromosome(self, chrom, chrom_df):
        # get the cur tiles for the chromosome
        cur = self.chrom[chrom]
        
        # gives each tiles' info an integer to map back to it
        tmp_to_info = {index: data[2] for index, data in enumerate(cur)}

        # converts the tiles' info to the integer
        cur = [i[:2] + [index] for index, i in enumerate(cur)]

        # convert the chrom_df to a list
        chrom_segs = chrom_df.apply(lambda x: [x.start_cm, x.end_cm, (x.id2, x.id1_haplotype, x.name)], axis=1).values.tolist()

        # add the new segments to the current tiles
        new_tile = segment_coverage(cur + chrom_segs, -self.gap)

        # iterate through the tiles and update the information is necessary
        for tile in new_tile:
            info = set()
            # the information is in tile[2] and is a set of (1) tuples and (2) index; must convert the index back to its info
            for i in tile[2]:
                info |= tmp_to_info.get(i, {i})
            
            tile[2] = info

        self.chrom[chrom] = new_tile

    # adds all segments from the segment_df df
    def create_tiles(self, segment_df):
        # add one chromosome at a time
        for chrom, chrom_df in segment_df.groupby("chromosome"):
            # add chrom segments
            self.add_chromosome(chrom, chrom_df)

    # returns the tile_df at the current state
    def return_tiles(self):

        tile_df = pd.DataFrame(columns=["coord", "segs1", "cluster1", "segs2", "cluster2"])

        for chrom, chrom_tile in self.chrom.items():
            # make a pandas df
            tmp = pd.DataFrame(chrom_tile, columns=["start_cm", "end_cm", "info"])

            # separate info into the two haplotypes' info
            tmp["info0"] = tmp["info"].apply(lambda x: [i for i in x if i[1]==0])
            tmp["info1"] = tmp["info"].apply(lambda x: [i for i in x if i[1]==1])
            
            # iterate through the haplotypes
            for hap in [0, 1]:
                tmp[f"cluster{hap+1}"] = tmp[f"info{hap}"].apply(lambda x: [i[0] for i in x])
                tmp[f"segs{hap+1}"] = tmp[f"info{hap}"].apply(lambda x: [i[2] for i in x])

            # add the coord column
            tmp["coord"] = tmp.apply(lambda x: (x.start_cm, x.end_cm), axis=1)

            # add the chromosome col
            tmp["chromosome"] = chrom

            tile_df = pd.concat([tile_df, tmp[["coord", "chromosome", "segs1", "cluster1", "segs2", "cluster2"]]])

        return tile_df.reset_index(drop=True)


'''
Takes a networkx graph object and performs spectral clustering on it.
'''
class SpectralCluster:

    def __init__(self, G, store_graph=False):

        ### Set all attributes to the default
        self.worked = False
        self.cluster1, self.cluster2 = nx.Graph(), nx.Graph()
        self.eigenval1, self.eigenval2 = np.inf, np.inf
        self.eigenvec1, self.eigenvec2 = np.array([]), np.array([]) 

        # empty graph
        if len(G.edges) < 1:
            return
        
        # too many nodes; downsample
        if len(G.nodes) > 10000:
            # keep the 10000 nodes belonging to the closest relatives
            tmp = pd.DataFrame(G.nodes, columns=["id2", "n", "k", "l"]).sort_values(by=["k", "l"]).tail(10000)
            
            largest = max(nx.connected_components(G.subgraph(tmp.apply(tuple, axis=1).values)), key=len)

            G = G.subgraph(largest)

        # create the signed Laplacian
        A = nx.adjacency_matrix(G).toarray()
        D = np.diag(np.sum(np.absolute(A),axis=1))
        # laplacian
        L = D - A

        # bool if signed
        signed = (A < 0).sum() > 0

        # something didn't work
        try:

            # eigendecomposition
            eigval, eigvector = np.linalg.eigh(L)

            # create eigenval/vector pair for sorting
            eigs = sorted([[val, vector] for val, vector in zip(eigval, eigvector.T)], key=lambda x: x[0])

            # cluster on 1st if signed, 2nd if not
            cluster_index = 0 if signed else 1
            
            # store the first two non-redundant
            self.eigenval1, self.eigenvec1 = eigs[cluster_index]
            self.eigenval2, self.eigenvec2 = eigs[cluster_index+1]

            self.worked = True

        # it did not work; quit running
        except:
            return

        # store the graph
        if store_graph:
            self.g = G

        # get the clusters
        cluster_bool = np.array([self.eigenvec1[n_index] > 0 for n_index, _ in enumerate(G.nodes)])
        nodes = list(G.nodes)
        self.cluster1 = G.subgraph([nodes[i] for i in np.where(cluster_bool)[0]])
        self.cluster2 = G.subgraph([nodes[i] for i in np.where(~cluster_bool)[0]])

'''Takes as input a tile_df and a t1 and t2 threshold (and a boolean whether the focal is male).
Returns a networkx graph object.'''
def graph_from_tiles(tile_df, t1, t2, male=False):

    # given two sub-tiles, returns an edge list of negative edges
    def negative_edges(cluster1, cluster2):
        return [[node1, node2, -1] for node1, node2 in it.product(cluster1, cluster2)]
    
    # given a cluster, returns an edge list of positive edges
    def positive_edges(cluster):
        return [[node1, node2, 1] for node1, node2 in it.combinations(cluster, r=2)]

    # create graph
    g = nx.Graph()

    # first add the negative between-subtile edges
    edges = list(it.chain(*tile_df.apply(lambda x: negative_edges(x.cluster1, x.cluster2), axis=1).values))
    # now add the positive between within-subtile edges
    edges += list(it.chain(*tile_df.apply(lambda x: positive_edges(x.cluster1) + positive_edges(x.cluster2), axis=1).values))

    # add edges between all nodes on the X-chromosome if the ind is male
    tmp = tile_df[tile_df["chromosome"].apply(lambda x: x==23 and male)]
    edges += positive_edges(it.chain(*tmp.apply(lambda x: x.cluster1 + x.cluster2, axis=1).values))

    # create a dataframe of the nodes
    nodes = it.chain(*tile_df.apply(lambda x: [[i, j, ("cluster1", x.name), x.chromosome==23 and male] for i,j in zip(x.cluster1, x.segs1)]\
                            + [[i, j, ("cluster2", x.name), x.chromosome==23 and male] for i,j in zip(x.cluster2, x.segs2)], axis=1).values)
    node_df = pd.DataFrame(nodes, columns=["tuple_id", "seg_index", "subtile", "maternal"])
    node_df[["id2", "n", "k", "l"]] = node_df["tuple_id"].apply(lambda x: list(x)).values.tolist()

    # only want to consider segments that are long enough or belong to close enough relatives
    tmp = node_df[node_df.apply(lambda x: x.l > t1 and x.k > t2, axis=1)]
    # add the final edges between segments of the same nodes that pass the edge thresholds
    edges += list(it.chain(*[positive_edges(df["tuple_id"].values) for _,df in tmp.groupby("id2")]))

    # add the edges
    g.add_weighted_edges_from(edges)

    # add info to the nodes for clustering
    attrs = node_df.apply(lambda x: {"seg_index": x.seg_index, "maternal": x.maternal, "subtile": x.subtile, "l": x.l}, axis=1).values
    nx.set_node_attributes(g, {node: attr for node, attr in zip(node_df["tuple_id"], attrs)})

    # largest subgraph
    try:
        largest_subgraph = sorted(list(nx.connected_components(g)), key=len)[-1]
    # no edge
    except IndexError:
        return nx.Graph()

    # return the largest subgraph
    return g.subgraph(largest_subgraph)

'''For a focal individual, stores and manages the clustering of their IBD segments'''
class FocalCluster:

    def __init__(self, focal, tile_df, male_bool):

        # gets the length of the tile
        tile_df["l"] = tile_df["coord"].apply(np.diff, axis=0).apply(sum)
        # true if the columns are ordered; if ordered will be ordered [maternal, paternal]
        tile_df["sex"] = False
        # haplotype columns
        tile_df["haplotype1"] = 0
        tile_df["haplotype2"] = 1

        self.focal = focal
        self.tile_df = tile_df
        self.male = male_bool
        self.clusterings = {}
        self.scheme = False


    def cluster(self, t1, t2):
        # get the graph from the tiles
        G = graph_from_tiles(self.tile_df, t1, t2, self.male)

        # spectral cluster
        clust = SpectralCluster(G)

        # clustering worked
        if clust.worked:    
            # get the sub-tiles
            tiles1 = pd.DataFrame({node[1]["subtile"] for node in clust.cluster1.nodes(data=True)}, columns=["parent1","tile"])
            tiles2 = pd.DataFrame({node[1]["subtile"] for node in clust.cluster2.nodes(data=True)}, columns=["parent2","tile"])

            # create the cluster df which has three columns
            # col1: tile number
            # col2: the subtile corresponding to the subtile of cluster1
            # col3: the subtile corresponding to the subtile of cluster2
            # two cols for the counts of relatives in each
            # drop dups bc these arise when there are multiple subtiles in the same parental cluster
            clust_df = tiles1.merge(tiles2, on="tile", how="outer").drop_duplicates("tile", keep=False)

            # remove any tiles where both clusters share the same subtile
            clust_df = clust_df[clust_df.apply(lambda x: x.parent1 != x.parent2, axis=1)]

            # the first parent should correspond to cluster1, this returns a bool whether to switch cluster1 and cluster2
            swap = ~clust_df.apply(lambda x: x.parent1 == "cluster1" or x.parent2 == "cluster2", axis=1)

            if swap.shape[0] == 0:
                clust.worked = False
            else:
                # add the tile index
                tmp = pd.DataFrame(swap.values, index=clust_df["tile"].values, columns=[(t1, t2)])

        if not clust.worked:
            self.tile_df[(t1, t2)] = np.nan
            clust.tot_cov = 0
            clust.x1, clust.x2 = 0, 0
            self.clusterings[(t1, t2)] = clust
            return [np.inf, np.nan, t1, t2]

        # merge with the existing tile_df
        self.tile_df = self.tile_df.merge(tmp, left_index=True, right_index=True, how="left", suffixes=["","_x"])

        # add the clustering information
        clust.tot_cov = self.tile_df[[(t1, t2), "l"]].dropna()["l"].sum()

        # add how much is shared on the X-chromosome (males only)
        clust.x1 = sum([node[1]["l"] for node in clust.cluster1.nodes(data=True) if node[1]["maternal"]])
        clust.x2 = sum([node[1]["l"] for node in clust.cluster2.nodes(data=True) if node[1]["maternal"]])

        # delete high mem attributes
        del clust.cluster1; del clust.cluster2; del clust.eigenvec1; del clust.eigenvec2

        # add to the dict of clusterings
        self.clusterings[(t1, t2)] = clust

        return [round(clust.eigenval1, 7), 1 / clust.tot_cov, t1, t2]

    def get_cluster_attr(self, t1, t2, attr):

        if (t1, t2) in self.clusterings:
            clust = self.clusterings[(t1, t2)]
            return getattr(clust, attr)
        
        else:
            return np.nan

    def return_tiles(self, t1, t2=None):
        # cols to return
        return_cols = ["chromosome", "coord", "cluster1", "cluster2", "segs1", "segs2", "haplotype1", "haplotype2", "sex"]

        # can return just the tile_df
        if t1 == None:
            return self.tile_df
        
        try:
            clust = self.clusterings[(t1, t2)]
        # clustering with the given thresholds des not exist
        except KeyError:
            return self.tile_df[self.tile_df.chromosome==-1][return_cols]

        # get the tile_df subset
        tiles = self.tile_df.dropna(axis=0, subset=[(t1, t2)]).copy()

        # swap the columns where necessary such that cluster1 corresponds to parent1
        swap = tiles[(t1, t2)]
        tiles = column_swapper(tiles, swap, [["cluster1", "cluster2"],
                                             ["segs1", "segs2"],
                                             ["haplotype1", "haplotype2"]])

        # get the total ibd of each cluster on the X-chromosome
        cluster1_X = clust.x1
        cluster2_X = clust.x2

        # there are no X-chromosome IBD segments or the focal is female or there is X IBD shared with both clusters
        if (cluster1_X + cluster2_X < 1) or (cluster1_X > 0 and cluster2_X > 0):
            # indicate that the columns are not ordered
            tiles["sex"] = False

        # there IS ibd on the X-chrom and the focal is male
        else:
            # indicate that the columns are ordered
            tiles["sex"] = True
            # want cluster1 to correspond to maternal relatives
            # cluster2 is the maternal cluster; swap all columns
            if cluster2_X > cluster1_X:
                tiles = tiles.rename(columns={"cluster1": "cluster2",
                                              "cluster2": "cluster1",
                                              "segs1": "segs2",
                                              "segs2": "segs1",
                                              "haplotype1": "haplotype2",
                                              "haplotype2": "haplotype1"})

        return tiles[return_cols]

    ''' cluster_scheme clusters on at least 10 t1, t2 combinations. It works from the most restrictive
    thresholds to the least restrictive and therefore we should see that the coverage increases and the
    smallest eigenvalue increases as we progress through the loop. '''
    def cluster_scheme(self, gain_t, t2s):

        # thresholds to iterate through
        t1s = [10, 5]

        # keep track of the clusterings
        clusterings = []

        # iterate through the clusterings
        for t1 in t1s:
            # init the prev clusters' coverage at 1 cM
            prev_cov = 1

            for t2 in t2s:
                # can't cluster
                if t1*2 > t2:
                    continue

                # get the current cluster and add to the clusterings
                cur_clust = self.cluster(t1, t2)
                clusterings.append(cur_clust)

                # get the coverage and the cur eigenval0
                cur_cov = 1 / cur_clust[1] if cur_clust[1] else 0
                cur_eig = cur_clust[0]

                # get the cM gain
                gain = (cur_cov - prev_cov) / prev_cov

                # stop the clustering
                if cur_eig > 0 and gain < gain_t:
                    print(f"Diminished tile gains for {self.focal} ({gain}) for t1={t1} and t2={t2}")
                    break

                prev_cov = max(cur_cov, 1)

        # sort the clusterings by eigenval0 and tot_cov
        clusterings.sort()

        # set the optimal
        for i, (_, _, t1, t2) in enumerate(clusterings):
            self.clusterings[(t1, t2)].optimal = i==0

        self.scheme = True

    # iterates through all clusterings and returns a simple summary of the clusterings
    def cluster_summary(self):

        if not self.scheme:
            return pd.DataFrame(columns=["focal", "t1", "t2", "tot_cov", "eigenval1", "optimal"])

        data = []
        for (t1, t2), clust in self.clusterings.items():
            # add the row
            data += [[self.focal, t1, t2, clust.tot_cov, clust.eigenval1, clust.optimal]]

        out_df = pd.DataFrame(data, columns=["focal", "t1", "t2", "tot_cov", "eigenval1", "optimal"])
        
        return out_df
    
    # iterate through the clusterings and return the optimal
    def get_optimal(self, tile_df=False):
        for (t1, t2), clust in self.clusterings.items():
            if clust.optimal:
                clust.t1 = t1; clust.t2 = t2
                return self.return_tiles(t1, t2) if tile_df else clust



'''Runs the clustering for a focal individual'''
def runClustering(focal, male_bool, segment_df, descendant_t, gain_t, t2s):

    # for a given segment, returns which clusters' tiles it overlaps with
    # could return {}, {cluster1}, {cluster2}, {cluster1, cluster2}
    def overlaps_tile(start, stop, hap, chrom, tile_df):
        # keeps track of which clusters (if any) the segment overlaps with
        overlaps = set()

        # iterate through each tile in the chromosome
        for _, row in tile_df[tile_df.chromosome==chrom].iterrows():
            # overlaps with the tile by 5+ cM
            if min(stop, row["coord"][1]) - max(start, row["coord"][0]) > 5:
                # which cluster the segment shares its haplotype with
                cluster = "cluster1" if row["haplotype1"] == hap else "cluster2"
                # add to the overlaps set
                overlaps |= {cluster}

        return list(overlaps)[0] if len(overlaps) == 1 else (None if len(overlaps)==0 else "both")

    print(f"Clustering focal ID {focal}")

    # init the tiles
    tiles = Tiles()

    # run the descendant finder
    if descendant_t > 0:

        # bool if the segment belongs to a putative descendant
        try:
            descendant = segment_df["id2"].apply(lambda x: x[2]) > descendant_t
        except:
            descendant = segment_df["id2"].apply(lambda x: x[2]) > descendant_t

        # split the segments into descendant/non-descendant dfs
        closeRel = segment_df.copy()[descendant]
        distantRel = segment_df.copy()[~descendant]

        # first create tiles with non-descendants
        tiles.create_tiles(distantRel)
        tile_df = tiles.return_tiles()

        # cluster using only the distant relatives' tiles
        clusters = FocalCluster(focal, tile_df, male_bool)
        clusters.cluster(10, 50)
        cluster_tiles = clusters.return_tiles(10, 50)

        if cluster_tiles.shape[0] == 0:
            clusters = FocalCluster(focal, tile_df, male_bool)
            clusters.cluster_scheme(gain_t, t2s)
            return clusters

        # add a col to each segment of the close relatives indicating which cluster they belong to
        closeRel["cluster"] = closeRel[["start_cm", "end_cm", "id1_haplotype", "chromosome"]].apply(lambda x: overlaps_tile(*x, cluster_tiles), axis=1)

        # add column that is just the relative identifier
        closeRel["id_name"] = closeRel["id2"].apply(lambda x: x[0])

        # add column that describes the descendant relationship if known
        closeRel["descendant"] = closeRel["id_name"].apply(lambda x: ("nibling" if "desc1" in x else "grandchild") if "desc" in x else "unknown")

        # iterate through the close relatives and keep track of the indexes to keep
        keep_segments = []
        for id2, df in closeRel.groupby("id_name"):

            desc_type = df["descendant"].values[0]

            tots = {c: (c_df["end_cm"] - c_df["start_cm"]).sum() for c, c_df in df.groupby("cluster")}

            # total IBD with each cluster
            clust1 = tots.get("cluster1", 0)
            clust2 = tots.get("cluster2", 0)

            # descendant
            if min(clust1, clust2) > 0.25*max(clust1, clust2):
                print(f"Found putative descendant for focal {focal} ({id2} shares {clust1} cM w/ cluster1 and {clust2} cM w/ cluster2). The relationship is a {desc_type} relationship. Total tile cov is {clusters.get_cluster_attr(10, 50, 'tot_cov')} cM")
            else:
                keep_segments += list(df.index)

                # false negative descendant
                if desc_type != "unknown":
                    print(f"Missed a descendant for focal {focal} ({id2} shares {clust1} cM w/ cluster1 and {clust2} cM w/ cluster2). This relative is a {desc_type} relationship. Total tile cov is {clusters.get_cluster_attr(10, 50, 'tot_cov')} cM")

        # now updates the tiles with the segments of the close relatives
        tiles.create_tiles(closeRel.loc[keep_segments])

    # Don't run the descendant-finder
    else:
        tiles.create_tiles(segment_df)

    # get the tile df
    tile_df = tiles.return_tiles()

    # now we can start clustering
    clusters = FocalCluster(focal, tile_df, male_bool)
    clusters.cluster_scheme(gain_t, t2s)

    # return clustering obj
    return clusters

class StoreResults:

    def __init__(self):

        self.focal_ids = []

    def add_focal(self, focal, cluster_obj):
        # add to list of focal ids
        self.focal_ids.append(focal)

        # add the clusering as an attribute
        setattr(self, focal, cluster_obj)

    # returns an iterator of all iids and their objs
    def iterator(self):
        return iter([(iid, getattr(self, iid)) for iid in self.focal_ids])

### returns the list of IBD files to read in
### if there are not batches, just returns the ibd file in a list
### otherwise returns a list of all possible combinations of batches
def relative_ibd_file_names(batches, relative_ibd_file):
    # no batches indicates, just return the ibd file
    if batches == ",,":
        return [relative_ibd_file]

    # first comma delim is the text delimiter for the batch
    text_delim = batches.split(",")[0]

    # cur batch is the batch of the current focal individuals
    cur_batch = int(batches.split(",")[1])

    # keeps track of which batch numbers to look at
    batch_list = []
    for b in batches.split(",")[2:]:
        if "-" in b:
            batch_list += [ i for i in range(int(b.split("-")[0]), int(b.split("-")[1]) + 1) ]
        else:
            batch_list.append(int(b))
    
    # writes the file names
    out_files = []
    for b1, b2 in it.product( [cur_batch], sorted(batch_list, key=lambda x: abs(x-cur_batch))):
        b1, b2 = sorted((b1, b2))
        batch1, batch2 = f"{text_delim}{b1}", f"{text_delim}{b2}"
        filename = relative_ibd_file.split(f"{text_delim}{cur_batch}")
        filename = filename[0] + batch1 + filename[1] + batch2 + filename[2]
        out_files.append(filename)
    return out_files

def runHaplotypeAnnotation(segments, duos, males, out):
    # init the hap df
    hap_df = pd.DataFrame()

    # have to run by chromosome
    for chrom, chrom_df in segments.groupby("chromosome"):
        # returns a dict of {focal: focal_df}
        hap_dfs = annotate_haplotypes(segments, duos, males)

        # have to concat the dfs and add the focal and chromosome cols
        for focal, focal_df in hap_dfs.items():
            focal_df["focal"] = focal
            focal_df["chromosome"] = chrom
            # concat to the existing df
            hap_df = pd.concat([hap_df, focal_df])

    # write out
    hap_df.reset_index(drop=True).to_feather(f"{out}_haplotypes.f")


if __name__ == "__main__":

    args = parse_args()
    sys.stdout.write("Parsed arguments.\n")

    ### read in batch information
    focal_file = pd.read_csv(args.focal_file, names = ["id1", "p1", "p2"], delim_whitespace=True, engine="python", dtype=str)
    sys.stdout.write("Read in the focal individual file.\n")

    # list of focal ID numbers
    args.focal_ids = focal_file["id1"].values.tolist()

    # list of duos
    args.duos = focal_file.dropna(subset=["p1"])[["id1","p1"]].apply(lambda x: tuple(x), axis=1).values.tolist()
    args.duos += focal_file.dropna(subset=["p2"])[["id1","p2"]].apply(lambda x: tuple(x), axis=1).values.tolist()

    # add the sexes
    args.sexes = {i[0]: i[1] for i in np.loadtxt(args.sex_file, dtype = str)} if args.sex_file != None else {}

    # get the list of relative ibd files
    args.relative_ibd_files = relative_ibd_file_names(args.batches, args.relative_ibd_file)

    # keep track of various times
    time_checkp = [time.time()]

    ### Load the IBD
    if args.ibd_pkl == "":

        IBD = IBDHandler(args)

        # each IBD file has all chromosomes in it
        if args.single_ibd_file:
            # load the IBD files
            IBD.load_ibd_from_list(args.relative_ibd_files, multiprocess=args.multiprocess)

        # separate chromosome for each
        else:
            for ibd_file in args.relative_ibd_files:
                # run a set of autosomal IBD files at a time
                IBD.load_ibd_from_list([ibd_file.replace("chr1", f"chr{chrom}") for chrom in range(1, 24)], args.multiprocess)

        if args.hap_annotation:
            # empty df to start
            segments = pd.DataFrame()

            # have to concat all the dfs but need to add the id1 col first
            for focal, focal_df in IBD.ibd.items():
                focal_df["id1"] = focal
                segments = pd.concat([segments, focal_df])

            # get list of male ids
            males = [id1 for id1, sex in args.sexes.items() if sex=="1"]

            runHaplotypeAnnotation(segments, args.duos, males, args.out)

            print("Haplotype annotation only. Exiting now...")
            sys.exit()

        # process the IBD segments
        IBD.process_all_ibd(args.multiprocess)

        time_checkp.append(time.time())
    
        out = open(f"{args.out}.ibd", "wb")
        pkl.dump(IBD, out)
        out.close()
        if args.ibd_only:
            print(f"IBD-only run. Exiting program now. Total run time: {round(time_checkp[1]-time_checkp[0], 3)}s")
            sys.exit()

    # ibd pkl file has been supplied
    else:
        infile = open(args.ibd_pkl, "rb")
        IBD = pkl.load(infile)
        infile.close()

        time_checkp.append(time.time())


    ### Run the clustering
    cluster_results = StoreResults()

    # load an existing results obj
    if args.continue_run:
        try:
            i = open(f"{args.out}_results.pkl", "rb")
            cluster_results = pkl.load(i)
            print(f"{len(cluster_results.focal_ids)} are already clustered. Skipping these samples.")
        except:
            pass

    # get focal IDs for focals that have IBD segments and are included in focal_ids (helpful for batching after the IBD obj has been created)
    focal_ids = [focal for focal, focal_df in IBD.ibd.items() if focal_df.shape[0] > 0 and focal in args.focal_ids and focal not in cluster_results.focal_ids]

    # run args.thread at a time
    threads = args.threads
    for index in range(0, len(focal_ids), threads):

        focal_iter = focal_ids[index:index+threads]

        with concurrent.futures.ProcessPoolExecutor() as executor:
            # run multiprocessing
            if args.multiprocess:

                results = [executor.submit(runClustering, focal, args.sexes.get(focal, "2")=="1", IBD.ibd[focal], args.descendant_k, args.coverage_gain, args.t2) for focal in focal_iter]
                # collect finished results
                clusters = [i.result() for i in concurrent.futures.as_completed(results)]
            # no multiprocessing; run one foal at a time
            else:
                clusters = [runClustering(focal, args.sexes.get(focal, "2")=="1", IBD.ibd[focal], args.descendant_k, args.coverage_gain, args.t2) for focal in focal_iter]

        # iterate through the results and add to the results file
        for focal_clusters in clusters:

            # only if clustering worked
            if not focal_clusters.scheme:
                continue

            # get the focal ID
            focal = focal_clusters.focal


            # get the ibd segments and delete to reduce memory usage
            focal_segments = IBD.ibd[focal]
            del IBD.ibd[focal]

            focal_clusters.segments = focal_segments
            del focal_segments

            cluster_results.add_focal(focal, focal_clusters)

        # write out the results
        o = open(f"{args.out}_results.pkl", "wb")
        pkl.dump(cluster_results, o)
        o.close()

    time_checkp.append(time.time())

    # print out times
    print(f"IBD analysis: {round(time_checkp[1]-time_checkp[0], 3)}s")
    print(f"Clustering: {round(time_checkp[2]-time_checkp[1], 3)}s")