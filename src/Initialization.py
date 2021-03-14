from src.SharedFunctions import *

def PrepareReference(
		Reference: str,
		PipelineConfigFile: str,
		Logger: logging.Logger = DefaultLogger("/dev/null")) -> None:
	
	MODULE_NAME = "PrepareReference"
	
	# Logging
	for line in [f"Reference: {Reference}"]: Logger.info(line)
	
	# Options
	Env = json.load(open(PipelineConfigFile, 'rt'))["GATKCondaEnv"]
	
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

def PrepareCapture(
		InputBED: str,
		Reference: str,
		OutputBED: str,
		Logger: logging.Logger) -> None:
	
	MODULE_NAME = "PrepareCapture"
	
	# Logging
	for line in [f"Input BED: {InputBED}", f"Reference: {Reference}", f"Output BED: {OutputBED}"]: Logger.info(line)
	
	with tempfile.TemporaryDirectory() as TempDir:
		
		# Options
		Faidx = Reference + ".fai"
		GenomeBED = os.path.join(TempDir, "genome.bed")
		FilteredBED = os.path.join(TempDir, "filtered.bed")
		
		# Processing
		PrepareGenomeBED(
			Reference = Reference,
			GenomeBED = GenomeBED,
			Logger = Logger)
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.Filter",
			Command = f"bedtools intersect -a \"{GenomeBED}\" -b \"{InputBED}\" > \"{FilteredBED}\"",
			Logger = Logger)
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.Sort",
			Command = f"bedtools sort -faidx \"{Faidx}\" -i \"{FilteredBED}\" | sed -e \'s/$/\\t\\./\' > \"{OutputBED}\"",
			Logger = Logger)

def Ucsc2Gff3(
		dbName: str,
		InputUCSC: str,
		OutputGFF3: str,
		Reference: str,
		Logger: logging.Logger) -> None:
	
	MODULE_NAME = "Ucsc2Gff3"
	
	# Logging
	for line in [f"Name: {dbName}", f"Input UCSC db: {InputUCSC}", f"Output GFF3: {OutputGFF3}"]: Logger.info(line)
	
	# Options
	AnchorCols = ["#chrom", "chromStart", "chromEnd"]
	DataOrder = ["#chrom", "sample", "type", "chromStart", "chromEnd", "score", "strand", "phase", "attributes"]
	
	# Prepare faidx
	Faidx = pandas.read_csv(Reference + ".fai", sep='\t', header=None).assign(Tag="##sequence-region", Start=1)[["Tag", 0, "Start", 1]]
	Chroms = {value: index for index, value in enumerate(Faidx[0].to_list())}
	
	# Load data
	Data = pandas.read_csv(InputUCSC, sep='\t')
	AttributeCols = [item for item in Data.columns.to_list() if item not in AnchorCols]
	
	# Filter & sort intervals by reference
	Data = Data[Data["#chrom"].apply(lambda x: x in Chroms.keys())]
	Data["Rank"] = Data["#chrom"].map(Chroms)
	Data.sort_values(["Rank", "chromStart"], inplace=True)
	
	# Processing
	Attributes = Data[AttributeCols].apply(lambda x: "ID=" + json.dumps(x.to_dict()).encode('utf-8').hex(), axis=1)
	Data.drop(columns=AttributeCols, inplace=True)
	Data["chromStart"] = Data["chromStart"].apply(lambda x: 1 if x == 0 else x)
	Data = Data.assign(sample=dbName, type="region", attributes=Attributes,	score=".", strand=".", phase=".")[DataOrder]
	
	# Save
	with open(OutputGFF3, 'wt') as O:
		O.write(
			"##gff-version 3\n" +
			Faidx.to_csv(sep=' ', index=False, header=False) + 
			Data.to_csv(sep='\t', index=False, header=False))
