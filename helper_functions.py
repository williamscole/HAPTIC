import numpy as np
import itertools as it
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

'''region dict is a dict of tuples (coordindates) mapping to a list of attributes,
i.e., a list of relatives who share IBD in some region
This function adds a new region (segment) to region dict
If it doesnt overlap with any of the existing regions, it simply creates a new one
If it does overlap, it creates new fragments of the overlap, appending relatives as
necessary to the newly split regions
Example input:
region_dict  = {(0, 10): [["A","B"]]}
new_region = [5, 10, ["C"]]
Returns: {(0, 5): [['A', 'B']], (5, 10): [['A', 'B'], ['C']]}
'''

def split_regions(region_dict, new_region):
    
    # returns the overlap of 2 regions (<= 0 if no overlap)
    def overlap(region1, region2):
        start1, end1 = region1
        start2, end2 = region2
        return min(end1,end2) - max(start1,start2)

    # out region will be returned; is a dict of regions mapping to members of region    
    out_region = dict()

    # overlapped keeps track of all the regions that overlap with the new region
    overlapped = {tuple(new_region[:2]):[new_region[2]]}

    # iterate through the existing regions
    for region in sorted(region_dict):

        # if overlap
        if overlap(region, new_region[:2]) > 0:

            # the regions completely overlap, just add the member and return region dict
            if tuple(region) == tuple(new_region[:2]):
                region_dict[region] += [new_region[2]]
                return region_dict

            # bc the region overlaps, add it to overlapped
            overlapped[region] = region_dict[region]

        # no overlap, but add the region to the out_region dict
        else:
            out_region[region] = region_dict[region]

    # all the segments in overlapped overlap, so each consecutive pairs of coordinates in sites should/could have different members
    sites = sorted(set(it.chain(*overlapped)))

    # iterate thru consecutive sites
    for start, stop in zip(sites, sites[1:]):

        # get the members of the regions that overlap the consecutive sites
        info = [j for i, j in overlapped.items() if overlap((start, stop), i) > 0]

        # unpack the membership
        out_region[(start,stop)] = sorted(it.chain(*info))

    return out_region

'''Same inputs as above but does something differently with them
If two regions overlap (or have a gap smaller than the gap param),
they are merged, provided they contain the same "info" (same relatives)
Example input:
region_dict = {(0, 20): ["A"]}
new_region = [21, 30, "A"]
gap = 1
Output: {(0, 30): ["A"]}
Warning: if our input was new_region = [0, 20, "B"], we'd lose the "A" information'''
def merge_regions(region_dict, new_region, gap = 0):

    # returns overlap of two regions (<= 0 if no overlap)
    def overlap(region1, region2):
        start1, end1 = region1
        start2, end2 = region2
        return min(end1,end2) - max(start1,start2)

    # see above
    out_region = {}
    overlapped = {tuple(new_region[:2]):new_region[2]}

    # iterate through the existing regions
    for region in sorted(region_dict):

        # True if the two regions are within some gap of each other or overlap
        overlap_bool = overlap(region, new_region[:2]) >= -gap

        # True if the regions belong to the same relative, same haplotype, etc.
        info_bool = region_dict[region] == new_region[2]

        # only overlapped if the two booleans are true
        if overlap_bool and info_bool:
            overlapped[region] = region_dict[region]

        # otherwise just add the region
        else:
            out_region[region] = region_dict[region]

    # the smallest and largest sites in overlapped are the start, stop sites of the total overlapped segment
    sites = sorted(it.chain(*overlapped))

    # get the info of all the sites
    info = [j for _,j in overlapped.items()]

    # add the overlapped sites to out region
    out_region[(sites[0],sites[-1])] = new_region[2]

    return out_region

def merge_regions(regions, new_region, gap = 0):

    # returns overlap of two regions (<= 0 if no overlap)
    def overlap(start1, end1, start2, end2):
        return min(end1,end2) - max(start1,start2)

    # see above
    out_region = list()
    overlapped = set(tuple(new_region[:2]))

    # iterate through the existing regions
    for region in sorted(regions, key = lambda x: (x[0], x[1])):

        # True if the two regions are within some gap of each other or overlap
        overlap_bool = overlap(*region[:2], *new_region[:2]) >= -gap

        # True if the regions belong to the same relative, same haplotype, etc.
        info_bool = region[2] == new_region[2]

        # only overlapped if the two booleans are true
        if overlap_bool and info_bool:
            overlapped |= set(region[:2])

        # otherwise just add the region
        else:
            out_region.append(region)

    # the smallest and largest sites in overlapped are the start, stop sites of the total overlapped segment
    sites = sorted(overlapped)

    # add the overlapped sites to out region
    out_region.append((sites[0], sites[-1], new_region[2]))

    return out_region


'''The goal of this function is to assign a PofO for the haplotypic
regions of the proband.
ibdfile is the pandas feather file with the PO segments
duos is a list of PO duos (child, parent)'''
def annotate_haplotypes(df, duos, males = []):

    # for duos, want to take the IBD segs shared with parent 1, switch the haplotype index, and rename the "parent"
    def duo_segments(duo_df):

        duo_df = duo_df[~duo_df.id1.isin(males)]

        # rename the missing parent as [id]_x
        duo_df["id2"] = duo_df["id1"].apply(lambda x: f"{x}_x")

        # switch the haplotype index
        duo_df["id1_haplotype"] = duo_df["id1_haplotype"].apply(lambda x: {0:1, 1:0}[x])

        return duo_df

    # calls either split_regions or merge_regions
    def regionalize(segments, func):
        regions = {}
        for seg in segments:
            regions = func(regions, seg)
        return regions
    
    # keep only PO pairs
    df["pair"] = df.apply(lambda x: (x.id1,x.id2),axis=1)
    df = df[df["pair"].isin(duos)].copy()

    ### if no F has a parent, continue
    if df.shape[0] == 0:
        return {}

    ### keep track of whos in a duo vs trio
    # trio maps id to boolean if they are in a trio
    trio = df.groupby("id1")["id2"].apply(lambda x: len(set(x))==2).to_dict()

    # creates boolean in df for trio or not
    df["trio"] = df["id1"].apply(lambda x: trio[x])

    # adds the segments of missing parents of duos
    df = pd.concat([df, duo_segments(df[~df.trio].copy())])

    # each focal maps to their hap annotation
    out_df = {}

    # iterate through each focal, haplotype pair
    for focal, segments in df.groupby(["id1", "id1_haplotype"]):

        # get segments and sort by the last coord of the segment
        segments = segments[["start_cm","end_cm","id2"]].values.tolist()
        segments.sort(key = lambda x: x[2])

        # first merge regions corresponding to the same parent
        merged_r = regionalize(segments, merge_regions)

        # merged_r is list of triple tuples: (start, stop, parent)
        # now merged between parents on the same haplotype
        split_r = regionalize(merged_r, split_regions)

        # create new df
        df = pd.DataFrame([[focal[0], focal[1], region[0], region[1], parents] for region, parents in split_r.items()],
                            columns=["F", "hap", "start", "stop", "parents"])

        # append df to other focal df
        out_df[focal[0]] = pd.concat([out_df.get(focal[0], pd.DataFrame()), df])

    return out_df

'''similar function to merge_regions except it doesnt require that the "info" of the two regions matches
Example input:
regions = [ [0, 10, "A"], [5, 10, "B"], [12, 20, "C"] ]
Output: [ [0, 10, ["A", "B"]], [12, 20, ["C"]]'''
def segment_coverage(regions, gap=0):
    #returns overlap
    def overlap(x1, y1, x2, y2):
        return min(y1,y2) - max(x1,x2)

    # each region has start, stop, set
    regions = [i[:2] + [set(i[2:3])] for i in regions]

    # to be returned
    out_regions = []

    # iterate through the regions
    for start, stop, d in regions:

        # for each iteration in regions, we iterate through out_regions to check for overlaps
        # overlapped coord keeps track of all overlapped coords between the region and out_regions
        # temp is the new out_region
        overlapped_coord, temp = [start, stop], []

        # iterate through out_regions
        for r1, r2, data in out_regions:

            # if overlapped/close enough gap, merge
            if overlap(start, stop, r1, r2) > -gap:
                overlapped_coord += [start, stop, r1, r2]

                # add the data
                d |= data

            # no overlap
            else:
                temp.append([r1,r2,data])

        # add the overlapped segments
        temp.append([min(overlapped_coord),max(overlapped_coord),d])

        # reassign out_regions
        out_regions = temp

    return out_regions

### for ukb with haplotype file
def parent_lab_dict(df, hap_df):

    def overlap(x1, y1, x2, y2):
        return min(y1, y2) - max(x1, x2)

    # returns the most "purple" node
    def min_purple(labs):
        min_purple = 0

        for _, d in labs.items():
            m = min(d[0], d[1])
            min_purple = m if m > min_purple else min_purple

        return min_purple

    def assign_parent(segments, regions, labs):
        
        regions = sorted(regions[["start","stop","parents"]].values.tolist())

        for index, row in segments.iterrows():

            id2, start, stop = row["id2"], row["start_cm"], row["end_cm"]

            p_kinship = {parent1: 0 , parent2: 0}

            # iterate through the parent haplotypes
            for r1, r2, parent in regions:

                # regions are sorted, stop loop if out of the window
                if r1 > stop:
                    break

                # compute the overlap between the segment and the region
                olap = overlap(start, stop, r1, r2)
                
                # if negative overlap, segment comes before the region, continue
                if olap <= 0:
                    continue

                # there is some overlap with the region; add cM to each parent of the region
                for p in parent:
                    p_kinship[p] += olap

            # take len of the segment
            l = stop - start

            # require that at least 80% of the segment is shared with a given parental hap
            p1 = p_kinship[parent1] / l > 0.8
            p2 = p_kinship[parent2] / l > 0.8

            # boolean pair to parent of origin of segment
            label = {(1, 1): 2,(1, 0): 0,(0, 1): 1,(0, 0): 3}[(p1, p2)]

            segment_labs[index] = label

            # add the segment
            if id2 not in labs:
                labs[id2] = {0: 0, 1: 0, 2: 0, 3: 0}
            labs[id2][label] += l

        return labs

    # given ibd sharing across whole genome, assigns the PofO
    def parent_label(d):
        pbool = (d[0] > 5, d[1] > 5)

        # if the amount of purple is <5 cM, otherwise it's purple
        pbool = pbool if d[2] < 5 else (1,1)

        return {(1, 1): 2, (1, 0): 0, (0, 1): 1, (0, 0): 3}[pbool]

    multi_chrom = [i for i, j in Counter(df.id2).items() if j > 1]
    
    parent1 = set(it.chain(*hap_df["parents"]))

    labs, segment_labs = {}, {}

    if len(parent1) == 2:

        parent1, parent2 = parent1

        for chrom, chrom_df in df.groupby("chromosome"):

            for hap, temp1 in chrom_df.groupby("id1_haplotype"):

                temp2 = hap_df[(hap_df["chromosome"] == chrom) & (hap_df["hap"] == hap)]

                labs = assign_parent(temp1, temp2, labs)


    return {id2: parent_label(d) for id2, d in labs.items()}, min_purple(labs), segment_labs

# for a given individual, returns purple node stats
# to be run in the assign_pofo function
def purple_stats(filename):

    if type(filename) == str:
        filename = [filename]

    ranges = [(10, 50), (50, 100), (100, 350), (350, 1260), (1260, 2100)]
    paired_data = {r: {cm: [0, 0] for cm in [5, 7.5, 10]} for r in ranges}
    agg_data = {cm: [0, 0] for cm in [5, 10, 15, 20, 50]}

    for n in filename:
        
        for focal, df in pd.read_feather(n).groupby("F"):
            print(focal)

            for index, row in df.iterrows():

                purple = [0, 0]

                k = 0

                for l, p in row["relative"]:

                    k += l
                    if p not in [0, 1]:
                        continue

                    purple[int(p)] += l

                for cm in [5, 10, 15, 20, 50]:
                    if purple[0] + purple[1] >= 2 * cm:
                        agg_data[cm][1] += 1

                    if purple[0] >= cm and purple[1] >= cm:
                        agg_data[cm][0] += 1

                for r, cm in it.product(ranges, [5, 7.5, 10]):

                    c = Counter(it.combinations([i[1] for i in row["relative"] if i[0] >= cm and r[0] <= k < r[1]], r = 2))
                    for i, j in c.items():
                                paired_data[r][cm][1] += j
                                paired_data[r][cm][0] += j if i[0] == i[1] else 0



    paired_data = pd.DataFrame([[r,c] + paired_data[r][c] for r, c in it.product(ranges, [5, 7.5, 10])], columns = ["range", "cM", "tp", "tot"])
    paired_data["p"] = paired_data["tp"] / paired_data["tot"]
    print(paired_data)

    for i, j in agg_data.items():
        print(i, j[0]/j[1])

    sns.catplot(x = "range", y = "p", hue = "cM", data = paired_data, kind = "bar")
    plt.savefig("test.png", dpi = 1000)



    
def get_roh_from_relative(df, overlap_cm):

    # sort the data frame
    sorted_df = df.sort_values(["chromosome", "id1", "id2", "id1_haplotype", "id2_haplotype", "start_cm"]).reset_index(drop=True)

    ### get the various booleans
    same_chrom = sorted_df.chromosome == sorted_df.shift(-1).chromosome
    same_id1 = sorted_df.id1 == sorted_df.shift(-1).id1
    same_id2 = sorted_df.id2 == sorted_df.shift(-1).id2
    diff_id1_hap = sorted_df.id1_haplotype != sorted_df.shift(-1).id1_haplotype
    same_id2_hap = sorted_df.id2_haplotype == sorted_df.shift(-1).id2_haplotype

    ### get the start and end of putative overlaps
    overlap_start = np.maximum(sorted_df.start_cm, sorted_df.shift(-1).start_cm)
    overlap_end = np.minimum(sorted_df.end_cm, sorted_df.shift(-1).end_cm)

    ### req all bools to be true and the overlap to be greater than overlap_cm
    roh_bool = (overlap_end - overlap_start > overlap_cm) & same_chrom & same_id1 & same_id2 & diff_id1_hap & same_id2_hap

    ### get the roh df
    roh_df = sorted_df[roh_bool].copy()

    ### set the start, end of the roh
    roh_df["start_cm"] = overlap_start[roh_bool]
    roh_df["end_cm"] = overlap_end[roh_bool]

    ### get self-roh and add it to the roh df
    self_ibd = sorted_df[sorted_df.id1 == sorted_df.id2].copy()
    roh_df = pd.concat([roh_df, self_ibd])

    ''' The issue we have now is that there are MANY segments corresponding to the same ROH.
    Here, we are merging these segments into a single ROH'''
    out_roh = []
    for chrom, chrom_df in roh_df.groupby("chromosome"):
        for id1, id1_df in chrom_df.groupby("id1"):
            # merge all ROH that overlap
            rohs = segment_coverage(id1_df[["start_cm", "end_cm"]].values.tolist(), 0)
            out_roh += [[id1, chrom] + i[:2] for i in rohs]
    
    out_roh = pd.DataFrame(out_roh, columns = ["id1", "chromosome", "start_cm", "end_cm"])
    out_roh["id2"] = out_roh["id1"]
    out_roh["id1_haplotype"] = 0
    out_roh["id2_haplotype"] = 1
    out_roh["roh"] = True

    # print(f"Found {out_roh.shape[0]} ROH.")

    return out_roh

'''Takes as input a pandas dataframe and a boolean Series that is True
if the columns need to swap. Swaps column pairs, which is a 2d list, where
each item is a pair of columns to switch'''
def column_swapper(df, bool, col_pairs):

    # create temp col
    df["temp_index"] = np.arange(df.shape[0])

    # store index
    index = df.index

    # create temp df
    tmp1 = df[bool]

    tmp2 = df[bool==False]

    for col1, col2 in col_pairs:
        tmp1 = tmp1.rename(columns={col1: col2, col2: col1})

    return pd.concat([tmp1, tmp2]).sort_values("temp_index").drop(["temp_index"], axis=1).set_index(index)
