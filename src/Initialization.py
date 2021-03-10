from src.SharedFunctions import *
from copy import deepcopy as dc

def PrepareReference(Reference, Logger, Env):
	
	MODULE_NAME = "PrepareReference"
	
	# Processing
	SimpleSubprocess(f"{MODULE_NAME}.SamtoolsIndex", f"samtools faidx \"{Reference}\"", Logger)
	SimpleSubprocess(f"{MODULE_NAME}.BWAIndex", f"bwa index \"{Reference}\"", Logger)
	SimpleSubprocess(f"{MODULE_NAME}.GATKIndex",  f"gatk CreateSequenceDictionary -R \"{Reference}\"", Logger, Env=Env) 

def Regenome(dbname, filename, dbtype, savename, Threads=cpu_count()):
	
	pandarallel.initialize(nb_workers=Threads, verbose=1)
	
	Tree = dict()
	dbtypes = ["ucsc"]
	
	if dbtype == dbtypes[0]:
		assert (type(dbname) == str) and (type(filename) == str), f"dbname and filename must be string type"
		chrom, start, end = "#chrom", "chromStart", "chromEnd"
		data = pandas.read_csv(filename, sep="\t")
		t1, t2 = data[[chrom, "chromStart"]].rename(columns={"chromStart": "pos"}), data[[chrom, "chromEnd"]].rename(columns={"chromEnd": "pos"})
		t1, t2 = t1.reset_index(), t2.reset_index()
		t1["type"], t2["type"] = True, False
		t = pandas.concat([t1, t2]).sort_values(by=[chrom, "pos"])
		stepper, now = [], {}
		for i, line in t.iterrows():
			if line["type"]: now[line["index"]] = True
			else: del now[line["index"]]
			stepper += [ dc(now) ]
		t["content"] = stepper
		newindex = [item for item in data.columns.to_list() if item not in [chrom, start, end]]
		t["content"] = t["content"].parallel_apply(lambda x: str([data.loc[item,:][newindex].to_dict() for item in x.keys()]))
		for c in list(set(data[chrom].to_list())):
			t0 = t[t[chrom] == c].sort_values(by="pos")
			Tree[c] = {"pos": t0["pos"].to_list(), "content": t0["content"].to_list()}
		SaveJSON(Tree, savename)
