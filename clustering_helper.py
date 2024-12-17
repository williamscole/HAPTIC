import numpy as np
import pandas as pd
import itertools as it

def focal_segments(focal, ibdfile = "ibdsegs/chr1.f"):
    def chrom_data(chrom):
        df = pd.read_feather(ibdfile.replace("1",str(chrom)))
        df1 = df[df["id1"]==focal]
        df2 = df[df["id2"]==focal]
        df2[["id1","id2"]] = df2[["id2","id1"]]
        df = df[df["id2"].apply(lambda x: "-" not in x)]
        return df
    df = chrom_data(1)
    for chrom in range(2,23):
        df = pd.concat([df, chrom_data(chrom)])
    df.reset_index(drop=True).to_feather("focal_dfs/ibdsegs_%s.f" % focal)

def focal_kinship(trio_list, batch_file = "focal_dfs/ukb_batch1.txt"):
    def pairs(batches, focal, pairs_to_keep):
        for b1, b2 in it.combinations_with_replacement(sorted(batches), r=2):
            print(b1,b2)
            batch1, batch2 = batches[b1], batches[b2]
            if b1 != b2:
                pairs = it.product(batch1,batch2)
            else:
                pairs = list(it.product(batch1,batch2)) + list(it.product(batch2,batch1))
            for pair in pairs:
                pairs_to_keep[(b1,b2)][pair] = pairs_to_keep[(b1,b2)].get(pair,[]) + [focal]
        return pairs_to_keep
    def add_kinship(row, b1, b2, pair):
        for focal in pairs_to_keep[(b1,b2)].get(pair, []):
            focal_kinship[focal].append(row.values.tolist())
    pairs_to_keep = {tuple(i):{} for i in it.combinations_with_replacement(np.arange(1,26),r=2)}
    focal_kinship = {}
    for focal,_,_ in trio_list:
        # df = pd.read_feather("../focal_dfs/ibdsegs_%s.f" % focal)
        relatives = set(df["id2"])
        batches = {i:set(np.loadtxt(batch_file.replace("1",str(i)),dtype=str)) & relatives for i in range(1,26)}
        pairs_to_keep = pairs(batches, focal, pairs_to_keep)
        focal_kinship[focal] = []
    for b1, b2 in it.combinations_with_replacement(np.arange(1,25)):
        kin = pd.read_feather("kinship/ALLchrom_batch%s_batch%s.kin" % (b1, b2))
        kin.apply(lambda x: add_kinship(x, b1, b2, tuple(x[0])), axis=1)

