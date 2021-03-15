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
