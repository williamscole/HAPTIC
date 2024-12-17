import pandas as pd
import pickle as pkl
import sys
import math
from relative_clustering import *
import warnings
import argparse
import numpy as np

# takes as input a vcf file; two matrices (both LxN) one is bool of switch or not, other is dist; list of sample ids (len of N)
# prints out the vcf
def open_switch_vcf(vcff, switch_m, tiled_m, region_m, sample_ids):

    # keep track of the index of the site
    site_index = 0

    # open the vcf and read line by line
    with open(vcff) as vcf:
        for line in vcf:

            # header line
            if "#" in line:

                # get the ids and rearrange the phase matrix
                if "#CHROM" in line:
                    # get the id and the associated index
                    id_list = [(i,index) for index, i in enumerate(line.split()[9:]) if i in sample_ids]
                    
                    # get the index of the ids in the phase matrix
                    idx = [sample_ids.index(i) for i,_ in id_list]

                    # rearrange the columns
                    switch_m = switch_m[:, idx]
                    tiled_m = tiled_m[:, idx]
                    region_m = region_m[:, idx]

                    # remove the ids in the vcf that are not phased and print the line
                    sample_idx = [index for _,index in id_list]
                    print("\t".join(line.split()[:9] + list(np.array(line.split()[9:])[sample_idx])))

                    continue

                # print all other header lines
                print(line.strip())

                # add the header for the phase type (TD)
                if "fileformat" in line:
                    print('##FORMAT=<ID=TD,Number=1,Type=Integer,Description="Describes how the site was phased. If -1, site is homozygous; if 0, site is phased in a tile; if >0, is the distance to the nearest tile">')
                    print('##FORMAT=<ID=RG,Number=1,Type=Integer,Description="The number of the region the site belongs to. If positive, part of a tiled region, else negative.">')
                continue

            # get the list of genotypes
            gts = [i[:3] for i in np.array(line.split()[9:])[sample_idx]]

            # create the new genotypes
            try:
                new_gts = [(f"{a2}{'|' if tiled else '/'}{a1}" if switch else f"{a1}{'|' if tiled else '/'}{a2}") + f":{str(tiled)}:{str(region)}" for (a1,_,a2), switch, tiled, region in zip(gts, switch_m[site_index].A1, tiled_m[site_index].A1, region_m[site_index].A1)]

                # replace the GT format field and print out the genotypes
                new_line = line.split()[:9]
                new_line[8] = "GT:TD:RG"
                print("\t".join(new_line + new_gts))
            except:
                print(site_index, switch_m.shape, tiled_m.shape)

            # next site index
            site_index += 1

def inmem_phase(vcff, switch_mat, id_list):
    import allel

    # open the VCF
    vcf = allel.read_vcf(vcff, fields="*", samples=id_list)

    # get the list of ids from the vcf
    vcf_samples = list(vcf["samples"])
    
    # get the indices of the vcf samples
    idx = [vcf_samples.index(i) for i in id_list]

    # m x n x 2 matrix
    gt_matrix = vcf["calldata/GT"][:,idx,:]

    # switch sites that need to be switched
    return np.where(switch_mat[:,:,np.newaxis], gt_matrix[:,:,[1,0]], gt_matrix)


def tile_to_switches(tile_df, map_df, chrom, inter_default, trim_tiles):
    # convert to str
    tile_df["chromosome"] = tile_df["chromosome"].apply(lambda x: 23 if x=="X" else int(x))

    # subset the tiles and add col for bool that's True if tiled
    chrom_df = tile_df[tile_df.chromosome==chrom].sort_values("coord").reset_index(drop=True)

    # no tiles; should return all sites as being no switch and -500
    if chrom_df.shape[0] == 0:
        return [True for _ in range(map_df.shape[0])], [-500 for _ in range(map_df.shape[0])], [-1 for _ in range(map_df.shape[0])]
    
    chrom_df["coord"] = chrom_df["coord"].apply(lambda x: (x[0]+trim_tiles, x[1]-trim_tiles))

    # add bool for tiled
    chrom_df["tiled"] = True

    # want cluster1 to correspond to haplotype 0; so wherever cluster1 is on haplotype1, we want ot switch
    chrom_df["switch"] = chrom_df["haplotype1"] == 1

    # add coord cols
    chrom_df["start"], chrom_df["end"] = zip(*chrom_df["coord"])

    # iterate through each consecutive pair of tiles
    for i in range(0, chrom_df.shape[0]-1):
        # overlap between the tile and the next
        if chrom_df.at[i, "end"] > chrom_df.at[i+1, "start"]:
            # get the start of the next, end of the cur
            start, end = chrom_df.at[i+1, "start"], chrom_df.at[i, "end"]
            # different switch status
            if chrom_df.at[i, "switch"] != chrom_df.at[i+1, "switch"]:
                # truncate the end of the first tile
                chrom_df.at[i, "end"] = start
            # truncate the start of the next tile
            chrom_df.at[i+1, "start"] = end

    chrom_df = chrom_df[chrom_df.apply(lambda x: x.end-x.start > 2.5, axis=1)].reset_index(drop=True)

    # get the inbetween regions; takes the phase switch of the previous tile
    tmp = pd.DataFrame({"start":chrom_df.end.shift(0), "end":chrom_df.start.shift(-1), "switch":chrom_df.switch}).fillna(1000)
    tmp = tmp[tmp.apply(lambda x: x.end > x.start, axis=1)]
    tmp["tiled"] = False

    # keep default
    if inter_default:
        tmp["switch"] = False

    # combine the regions
    region = [[-500, chrom_df.at[0, "start"], chrom_df.at[0, "switch"], False]]
    region += chrom_df[["start", "end", "switch", "tiled"]].values.tolist()
    region += tmp[["start", "end", "switch", "tiled"]].values.tolist()
    region.sort(key=lambda x: x[0])

    # iterate through the regions and subset the map file
    dfs = []
    region_n = 1
    for start, end, switch, tiled in region:
        tmp = map_df[map_df["CM"].apply(lambda x: start <= x < end)].copy()
        # no sites
        if tmp.shape[0] == 0:
            continue
        tmp["switch"] = switch
        # get the distance to the nearest tile start/end
        tmp["tiled"] = tmp["CM"].apply(lambda x: str(math.ceil(min(end-x, x-start))*(1 if tiled else -1)))
        tmp["region"] = region_n if tiled else -region_n
        region_n += 1
        dfs.append(tmp)

    # concat all the dataframes
    new_map = pd.concat(dfs)

    # return first the switch bool and then the distances
    return new_map["switch"].values.tolist(), new_map["tiled"].values.tolist(), new_map["region"].values.tolist()


### Assesses the whole-chromosome phase quality; requires that all unphased sites in the trio vcf are set to missing
def phase_quality(vcff, trio_vcff, output, tile_dist, **kwargs):
    def compute_phase_error(mat, samples, region):
        # get the number of het sites considered
        n_het = (~np.isnan(mat)).sum(axis=0)

        # get the number of errors
        n_err = np.nansum(mat, axis=0)

        # get the PER
        switch = np.divide(n_err, n_het, out=np.zeros_like(n_err), where=(n_het != 0))

        tmp = pd.DataFrame({"INDV": samples, "N_COMMON_PHASED_HET": n_het, "N_SWITCH": n_err, "SWITCH": switch})
        
        tmp["region"] = region

        return tmp[tmp.N_COMMON_PHASED_HET>0]
    
    def regional_phase_error(mat, region_mat, samples):

        # samples that have no tiles on the chromosome and then set the genotype to nan
        dfs = [compute_phase_error(np.where(region_mat[-1,:]==-1, mat, np.nan), samples, "chromosome")]
        mat = np.where(region_mat[-1,:] != -1, mat, np.nan)

        # samples that start with a gap
        dfs.append(compute_phase_error(np.where(region_mat==-1, mat, np.nan), samples, "start"))
        mat = np.where(region_mat != -1, mat, np.nan)

        # samples that end with a gap
        dfs.append(compute_phase_error(np.where(region_mat==region_mat[-1,:], mat, np.nan), samples, "end"))
        mat = np.where(region_mat != region_mat[-1,:], mat, np.nan)

        # get the max number of tiles
        max_n_tiles = np.abs(region_mat[-1,:]).max()

        for n in range(1, max_n_tiles+1):
            for region, n in zip(["tiled", "gap"], [n, -n]):
                dfs.append(compute_phase_error(np.where(region_mat==n, mat, np.nan), samples, region))

        return pd.concat(dfs).sort_values("INDV").reset_index(drop=True)

    import allel

    # load the vcf of the stat phased/cluster phased
    if type(vcff) == str:
        vcf = allel.read_vcf(vcff, fields=["samples", "calldata/GT"])
    # the data have been phased inmem
    else:
        vcf = vcff
    samples = vcf["samples"]

    # load the vcf of the trio-phased vcf
    trio = allel.read_vcf(trio_vcff, fields=["samples", "calldata/GT", "variants/POS", "variants/CHROM"], samples=vcf["samples"])
    trio_samples = list(trio["samples"])

    # check to see that all samples have a corresponding trio sample
    no_trio_sample = [i for i in samples if i not in trio_samples]
    if len(no_trio_sample) > 0:
        m,_,_ = trio["calldata/GT"].shape
        nan_cols = np.full((m, len(no_trio_sample), 2), np.nan)
        trio["calldata/GT"] = np.concatenate([trio["calldata/GT"], nan_cols], axis=1)
        trio_samples += no_trio_sample

    # the index of each individual in the trio samples data
    idx = [trio_samples.index(i) for i in samples]

    # this is a m x n matrix with reordered columns
    trio_matrix = trio["calldata/GT"][:,idx,0]

    # this is a m x n x 2 matrix
    gt_matrix = vcf["calldata/GT"]

    # get all the het. sites
    het_sites = gt_matrix[:,:,0]!=gt_matrix[:,:,1]

    # get all the sites that are otherwise missing
    nonmissing = (gt_matrix[:,:,0]>-1) & (gt_matrix[:,:,1]>-1) & (trio_matrix>-1)

    # if either homozyg or missing, set the site to NaN
    hap0_matrix = np.where(het_sites & nonmissing, gt_matrix[:,:,0], np.nan)

    d1, d2 = tile_dist
    if d1 + d2 != 0:
        td = vcf["calldata/TD"]
        hap0_matrix = np.where((td >= d1) & (td <= d2), hap0_matrix, np.nan)

    # get the difference in genotype; for sites not considered (which are nan) this will return nan
    mat = np.abs(trio_matrix - hap0_matrix)

    if kwargs["region_mat"].shape[0] != 0:
        out_df = regional_phase_error(mat, kwargs["region_mat"], samples)

    else:
        out_df = compute_phase_error(mat, samples, "none")

    out_df.to_csv(output, index=False, sep="\t")




### Takes all chromosomes and concats
def concat_chroms(diff_file):

    diff_df = pd.read_csv(diff_file, delim_whitespace=True, header=None, names=["INDV","N_COMMON_PHASED_HET", "N_SWITCH", "SWITCH", "region"])

    print("INDV\tN_COMMON_PHASED_HET\tN_SWITCH\tSWITCH\tregion")

    for indv, indv_df in diff_df.groupby("INDV"):

        # get the number of sites
        n_sites = indv_df["N_COMMON_PHASED_HET"].sum()
        n_switch = indv_df["N_SWITCH"].sum()
        switch = n_switch / n_sites if n_sites > 0 else 0

        # take the complement phase error
        indv_df["N_SWITCH"] = indv_df.apply(lambda x: x.N_COMMON_PHASED_HET-x.N_SWITCH if switch > 0.5 else x.N_SWITCH, axis=1)

        # print switch error by region
        for region, region_df in indv_df.groupby("region"):
            for _, row in region_df.iterrows():
                print(f"{indv}\t{row['N_COMMON_PHASED_HET']}\t{row['N_SWITCH']}\t{row['N_SWITCH']/row['N_COMMON_PHASED_HET']}\t{region}")


def debug_args(args):
    args.vcf = f"/gpfs/data/sramacha/ukbiobank_jun17/cwilli50/relative-clustering/debug_phaser/chr{args.chr}_statPhased.vcf" if args.vcf == None else args.vcf
    args.trio = f"/gpfs/data/sramacha/ukbiobank_jun17/cwilli50/relative-clustering/debug_phaser/chr{args.chr}_trioPhased.vcf.gz" if args.trio == None else args.trio
    args.map = f"/gpfs/data/sramacha/ukbiobank_jun17/cwilli50/relative-clustering/debug_phaser/chr{args.chr}.map" if args.map == None else args.map
    args.results = f"/gpfs/data/sramacha/ukbiobank_jun17/cwilli50/relative-clustering/debug_phaser/results.pkl" if args.results == None else args.results
    args.output = f"/gpfs/data/sramacha/ukbiobank_jun17/cwilli50/relative-clustering/debug_phaser/tmp_chr{args.chr}.txt" if args.output == None else args.output
    args.add_underscore = True
    return args


def parse_args():
    """Parses the command line arguments and returns the ArgumentParser.parse_args() object."""
    parser = argparse.ArgumentParser()

    ### Basic arguments
    parser.add_argument("-add_underscore", action="store_true")
    parser.add_argument("-debug", action="store_true")
    parser.add_argument("-chr", type=str)

    ### Phase the VCF based on the tile
    parser.add_argument("-phase", help="Phase. Must use the -vcf flag, -results flag, and -map flag.", action="store_true")
    parser.add_argument("-vcf", help="vcf file for input.", type=str)
    parser.add_argument("-map", help="Path to .map file for chromosome 1", type=str)
    parser.add_argument("-results", help="Phasing results .pkl file", type=str)
    parser.add_argument("-output", help="Full path + file name for the output file")

    parser.add_argument("-write_phase", help="Must use with -phase if you want to print out the VCF.", action="store_true")
    parser.add_argument("-thresholds", nargs='+', help="The t1 and t2 thresholds (in that order).", default=[0, 0], type=int)
    parser.add_argument("-inter_default", action="store_true", help="If True, does not attempt to fix the phase.")
    parser.add_argument("-trim_tiles", help="In cM, the amount to trim off the ends of the tiles.", default=0.0, type=float)


    ### Assess the PER rate of the phased dataset
    parser.add_argument("-assess", help="Assess the phase quality. Must use the -trio flag and the -vcf flag.", action="store_true")
    parser.add_argument("-trio", help="Trio-phased file. Unphased sites must be set to missing.", type=str)
    parser.add_argument("-regional", action="store_true", help="Prints PER for each phased/gapped region of the genome.")
    # -vcf flag must be used here as well

    parser.add_argument("-tile_dist", nargs='+', help="The thresholds (inclusive) of distance to nearest tile (negative if not in a tile; positive if in a tile).", default=[-500, 500], type=float)

    ### Concat a bunch of .diff files
    parser.add_argument("-concat", help="Concat the files.", action="store_true")
    parser.add_argument("-diff", help=".diff file containing all chromosomes.")

    args = parser.parse_args()
    
    return args

if __name__ == "__main__":

    args = parse_args()

    if args.debug:
        args = debug_args(args)

    # phase
    if args.phase:

        # load results
        i = open(args.results, "rb")
        results = pkl.load(i)

        # load map file
        map_df = pd.read_csv(args.map, delim_whitespace=True, header=None, names=["CHROM", "ID", "CM", "POS"])

        # get the tiles and the new map dfs
        map_dfs, focal_ids = [], []
        for focal, r in results.iterator():

            # get the tile_df; either the optimal or from the specified thresholds
            tile_df = r.get_optimal(tile_df=True) if sum(args.thresholds)==0 else r.return_tiles(*args.thresholds)

            # add the phase/tile information
            map_dfs.append(tile_to_switches(tile_df, map_df, int(args.chr) if args.chr != "X" else 23, args.inter_default, args.trim_tiles))

            # add the focal
            focal_ids.append(f"{focal}_{focal}" if args.add_underscore else focal)

        # get the matrices
        switch_m = np.matrix([i for i,_,_ in map_dfs]).T
        tiled_m = np.matrix([j for _,j,_ in map_dfs]).T
        region_m = np.matrix([k for _,_,k in map_dfs]).T

        if args.write_phase:
            # phase
            open_switch_vcf(args.vcf, switch_m, tiled_m, region_m, focal_ids)
        else:
            gt_matrix = inmem_phase(args.vcf, switch_m, focal_ids)

    # cannot write_phase and assess in the same step; rerun without the -phase flag
    if args.assess and not args.write_phase:

        if "gt_matrix" in globals():
            phase_quality({"calldata/GT": gt_matrix, "samples": focal_ids, "calldata/TD": tiled_m.astype(int)},
                          args.trio, args.output, args.tile_dist, region_mat=region_m if args.regional else np.zeros((0, 0)))
            
        else:
            phase_quality(args.vcf, args.trio, args.output, args.tile_dist, region_mat=region_m if args.regional else np.zeros((0, 0)))

    elif args.concat:

        concat_chroms(args.diff)

