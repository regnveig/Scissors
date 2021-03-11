from src.SharedFunctions import *

def PrepareReference(
	Reference: str,
	Logger: logging.Logger,
	Env: str) -> None:
	MODULE_NAME = "PrepareReference"
	# Processing
	SimpleSubprocess(
		Name = f"{MODULE_NAME}.SamtoolsIndex",
		Command = f"samtools faidx \"{Reference}\"",
		Logger = Logger)
	SimpleSubprocess(
		Name = f"{MODULE_NAME}.BWAIndex",
		Command = f"bwa index \"{Reference}\"",
		Logger = Logger)
	SimpleSubprocess(
		Name = f"{MODULE_NAME}.GATKIndex",
		Command = f"gatk CreateSequenceDictionary -R \"{Reference}\"",
		Logger = Logger,
		Env = Env)
	SimpleSubprocess(
		Name = f"{MODULE_NAME}.GenomeBED",
		Command = "awk 'BEGIN {FS=\"\\t\"}; {print $1 FS \"0\" FS $2}' \"" + Reference + ".fai\" > \"" + Reference + ".bed\"",
		Logger = Logger)

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
