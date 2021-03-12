from src.SharedFunctions import *

## ------======| PIPELINE STAGES |======------

## Trimming & Convertation Reads

def Solid2Illumina(InputFQ: str,
				   OutputFQ: str,
				   Logger: logging.Logger) -> None:
	MODULE_NAME = "Solid2Illumina"
	# Logging
	for line in [
		f"Input FASTQ: {InputFQ}",
		f"Output FASTQ: {OutputFQ}"]: Logger.info(line)
	# Processing
	SimpleSubprocess(
		Name = f"{MODULE_NAME}.Convert",
		Command = f"cutadapt -c --format=sra-fastq --bwa --action=none -o \"{OutputFQ}\" \"{InputFQ}\"",
		Logger = Logger)
	
def Cutadapt(InputR1: str,
			 InputR2: str,
			 OutputR1: Union[str, None],
			 OutputR2: Union[str, None],
			 Adapter: dict,
			 ReportTXT: str,
			 Threads: int,
			 Logger: logging.Logger) -> None:
	MODULE_NAME = "Cutadapt"
	# Logging
	for line in ([f"Mode: Single-end", f"Input FASTQ: {InputR1}", f"Output FASTQ: {OutputR1}"] if InputR2 is None else [f"Mode: Paired-end", f"Input FASTQ [R1]: {InputR1}", f"Input FASTQ [R2]: {InputR2}", f"Output FASTQ [R1]: {OutputR1}", f"Output FASTQ [R2]: {OutputR2}"]) + [f"Adapter: {Adapter['Name']}"]: Logger.info(line)
	# Processing
	SimpleSubprocess(
		Name = f"{MODULE_NAME}.Trim",
		Command = f"cutadapt -j {str(Threads)} -e 0.2 -m 8 -a {Adapter['R1']} -o \"{OutputR1}\" \"{InputR1}\" > \"{ReportTXT}\"" if InputR2 is None else f"cutadapt -j {str(Threads)} -e 0.2 -m 8 -a {Adapter['R1']} -A {Adapter['R2']} -o \"{OutputR1}\" -p \"{OutputR2}\" \"{InputR1}\" \"{InputR2}\" > \"{ReportTXT}\"",
		Logger = Logger)

def FastQC(InputFastQ: str,
		   OutputHTML: str,
		   Logger: logging.Logger,
		   Size: int = 0,
		   Threads: int = cpu_count()) -> None:
	MODULE_NAME = "FastQC"
	# Logging
	for line in [f"FastQ: {InputFastQ}", f"Report: {OutputHTML}"]: Logger.info(line)
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		AnalyzeFilename = InputFastQ
		if Size != 0:
			Logger.info(f"Subsample: {str(Size)}")
			SampleFilename = os.path.join(TempDir, "sample.fastq.gz")
			SimpleSubprocess(
				Name = f"{MODULE_NAME}.Sampling",
				Command = f"zcat -q \"{InputFastQ}\" | head -{str(Size * 4)} | gzip -c > \"{SampleFilename}\"",
				Logger = Logger)
			AnalyzeFilename = SampleFilename
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.Analysis",
			Command = f"fastqc -o \"{TempDir}\" -t {str(Threads)} \"{AnalyzeFilename}\"",
			Logger = Logger)
		HTMLTemp = glob(os.path.join(TempDir, "*.html"))
		if len(HTMLTemp) != 1:
			ErrorMessage = f"Error processing file '{InputFastQ}'"
			Logger.error(ErrorMessage)
			raise RuntimeError(ErrorMessage)
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.Move",
			Command = f"cp \"{HTMLTemp[0]}\" \"{OutputHTML}\"",
			Logger = Logger)

## Alignment & Merge BAMs

def BWA(InputR1: str,
		InputR2: Union[str, None],
		Reference: str,
		RGHeader: str,
		OutputBAM: str,
		Logger: logging.Logger,
		Env: str,
		Threads: int = cpu_count()) -> None:
	MODULE_NAME = "BWA"
	# Logging
	for line in ([f"Reads: {InputR1}"] if InputR2 is None else [f"R1: {InputR1}", f"R2: {InputR2}"]) + [f"Output BAM: {OutputBAM}", f"Reference: {Reference}", f"RG Header: {RGHeader}"]: Logger.info(line)
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.AlignAndSort",
			Command = f"bwa mem -R \"{RGHeader}\" -t {str(Threads)} -v 1 \"{Reference}\" \"{InputR1}\" | gatk SortSam --VERBOSITY ERROR --TMP_DIR \"{TempDir}\" -SO queryname -I \"/dev/stdin\" -O \"{OutputBAM}\"" if (InputR2 is None) else f"bwa mem -R \"{RGHeader}\" -t {str(Threads)} -v 1 \"{Reference}\" \"{InputR1}\" \"{InputR2}\" | gatk SortSam --VERBOSITY ERROR --TMP_DIR \"{TempDir}\" -SO queryname -I \"/dev/stdin\" -O \"{OutputBAM}\"",
			Logger = Logger,
			CheckPipefail = True,
			Env = Env)

def MergeBAMs(BAMs: list,
			  OutputBAM: str,
			  SortOrder: str,
			  Logger: logging.Logger,
			  Env: str) -> None:
	MODULE_NAME = "MergeBAMs"
	# Processing
	SimpleSubprocess(
		Name = f"{MODULE_NAME}.Merge",
		Command = f"gatk MergeSamFiles --USE_THREADING true -SO {SortOrder} {MultipleTags(Tag='-I', List=BAMs)} -O \"{OutputBAM}\"",
		Logger = Logger,
		Env = Env)

## Mark Duplicates

def MarkDuplicates(InputBAM: str,
				   OutputBAM: str,
				   MetricsTXT: str,
				   Logger: logging.Logger,
				   Env: str) -> None:
	MODULE_NAME = "MarkDuplicates"
	# Logging
	for line in [f"Input: {InputBAM}", f"Output: {OutputBAM}", f"Metrics: {MetricsTXT}"]: Logger.info(line)
	# Options
	JavaOptions = f"-XX:+UseParallelGC -XX:ParallelGCThreads=2" # -Xmx12G removed, can't see the sense of memory limits if the process is single-threaded.
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.RemoveAndSort",
			Command = f"gatk --java-options \"{JavaOptions}\" MarkDuplicates --REMOVE_DUPLICATES true --VERBOSITY ERROR --ASSUME_SORT_ORDER queryname --TMP_DIR \"{TempDir}\" -M \"{MetricsTXT}\" -I \"{InputBAM}\" -O \"/dev/stdout\" | gatk SortSam --VERBOSITY ERROR --TMP_DIR \"{TempDir}\" -SO coordinate -I \"/dev/stdin\" -O \"{OutputBAM}\"",
			Logger = Logger,
			CheckPipefail = True,
			Env = Env)

## BQSR

def ContigBaseRecalibration(Contig: str,
							InputBAM: str,
							TempDir: str,
							dbSNP: str,
							Reference: str,
							GATK_ConfigFile: str,
							Logger: logging.Logger,
							Env: str) -> str:
	MODULE_NAME = "ContigBaseRecalibration"
	# Options
	FiltersComparison = ["MappedReadFilter", "MappingQualityAvailableReadFilter", "MappingQualityNotZeroReadFilter", "NotDuplicateReadFilter", "NotSecondaryAlignmentReadFilter", "PassesVendorQualityCheckReadFilter"] 
	# This is a bug above. BaseRecalibrator and ApplyBQSR filter reads differently,
	# so ApplyBQSR f***s up every time it can't find RG.
	BaseRecalibrator_JavaOptions = f"-Xmx3G -XX:+UseParallelGC -XX:ParallelGCThreads=2"
	ApplyBQSR_JavaOptions = f"-Xmx3G"
	BQSRTable = os.path.join(TempDir, f"bqsr_table_{Contig}.tsv")
	OutputBAM = os.path.join(TempDir, f"output_{Contig}.bam")
	# Processing
	SimpleSubprocess(
		Name = f"{MODULE_NAME}.MakeTable[{Contig}]",
		Command = f"gatk --java-options \"{BaseRecalibrator_JavaOptions}\" BaseRecalibrator --gatk-config-file \"{GATK_ConfigFile}\" --tmp-dir \"{TempDir}\" -L {Contig} -I \"{InputBAM}\" --known-sites \"{dbSNP}\" -O \"{BQSRTable}\" -R \"{Reference}\"",
		Logger = Logger,
		Env = Env)
	SimpleSubprocess(
		Name = f"{MODULE_NAME}.Apply[{Contig}]",
		Command = f"gatk --java-options \"{ApplyBQSR_JavaOptions}\" ApplyBQSR --gatk-config-file \"{GATK_ConfigFile}\" {MultipleTags(Tag='-RF', List=FiltersComparison, Quoted=False)} --tmp-dir \"{TempDir}\" -OBI false -L {Contig} -bqsr \"{BQSRTable}\" -I \"{InputBAM}\" -O \"{OutputBAM}\"",
		Logger = Logger,
		Env = Env)
	# Return
	return OutputBAM

def BaseRecalibration(InputBAM: str,
					  OutputBAM: str,
					  dbSNP: str,
					  Reference: str,
					  GATK_ConfigFile: str,
					  Logger: logging.Logger,
					  Env: str,
					  Threads: int = cpu_count()) -> None:
	MODULE_NAME = "BaseRecalibration"
	# Logging
	for line in [f"Input: {InputBAM}", f"Output: {OutputBAM}", f"Known sites: {dbSNP}", f"Reference: {Reference}"]: Logger.info(line)
	# Options
	Contigs = [item['SN'] for item in GetContigs(FileBAM=InputBAM)]
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.PreIndex",
			Command = f"gatk BuildBamIndex -I \"{InputBAM}\"",
			Logger = Logger,
			Env = Env)
		with Threading("ContigBaseRecalibration", Logger, Threads) as pool:
			Shards = pool.map(functools.partial(
				ContigBaseRecalibration,
				InputBAM = InputBAM,
				TempDir = TempDir,
				dbSNP = dbSNP,
				Reference = Reference,
				GATK_ConfigFile = GATK_ConfigFile,
				Logger = Logger,
				Env = Env), Contigs)
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.Merge",
			Command = f"gatk MergeSamFiles --USE_THREADING true -SO coordinate {MultipleTags(Tag='-I', List=Shards)} -O \"{OutputBAM}\"",
			Logger = Logger,
			Env = Env)
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.PostIndex",
			Command = f"gatk BuildBamIndex -I \"{OutputBAM}\"",
			Logger = Logger,
			Env = Env)

## Haplotype Calling

def ContigHaplotypeCalling(Contig: str,
						   InputBAM: str,
						   TempDir: str,
						   Reference: str,
						   Logger: logging.Logger,
						   Env: str) -> str:
	MODULE_NAME = "ContigHaplotypeCalling"
	# Options
	JavaOptions = f"-Xmx3G"
	OutputVCF = os.path.join(TempDir, f"output_{Contig}.vcf")
	# Processing
	SimpleSubprocess(
		Name = f"{MODULE_NAME}.Calling[{Contig}]",
		Command = f"gatk --java-options \"{JavaOptions}\" HaplotypeCaller --native-pair-hmm-threads 2 -OVI false --dont-use-soft-clipped-bases true -L {Contig} -I \"{InputBAM}\" -O \"{OutputVCF}\" -R \"{Reference}\"",
		Logger = Logger,
		Env = Env)
	# Return
	return OutputVCF

def HaplotypeCalling(InputBAM: str,
					 OutputVCF: str,
					 Reference: str,
					 Logger: logging.Logger,
					 Env: str,
					 Threads: int = cpu_count()) -> None:
	MODULE_NAME = "HaplotypeCalling"
	# Logging
	for line in [f"Input BAM: {InputBAM}", f"Output VCF: {OutputVCF}", f"Reference: {Reference}"]: Logger.info(line)
	# Options
	Contigs = [item['SN'] for item in GetContigs(FileBAM=InputBAM)]
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		with Threading("ContigHaplotypeCalling", Logger, Threads) as pool:
			Shards = pool.map(functools.partial(
				ContigHaplotypeCalling,
				InputBAM = InputBAM,
				TempDir = TempDir,
				Reference = Reference,
				Logger = Logger,
				Env = Env), Contigs)
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.Merge",
			Command = f"gatk MergeVcfs {MultipleTags('-I', Shards)} -O \"{OutputVCF}\"",
			Logger = Logger,
			Env = Env)

## Statistics

def CoverageStats(Name: str,
				  FinalBAM: str,
				  StatsTXT: str,
				  CaptureBED: str,
				  Reference: str,
				  Logger: logging.Logger) -> None:
	MODULE_NAME = "CoverageStats"
	# Logging
	for line in [f"BAM File: {FinalBAM}", f"Reference: {Reference}", f"Capture BED: {CaptureBED}"]: Logger.info(line)
	with tempfile.TemporaryDirectory() as TempDir:
		# Options
		RefIndex = f"{Reference}.fai"
		NotCaptureBED = os.path.join(TempDir, "not_capture.bed")
		CaptureTemp = os.path.join(TempDir, "capture.csv")
		NotCaptureTemp = os.path.join(TempDir, "not_capture.csv")
		GenomeBED = os.path.join(TempDir, "genome.bed")
		# Coverage
		PrepareGenomeBED(
			Reference = Reference,
			GenomeBED = GenomeBED,
			Logger = Logger)
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.CreateNotCaptureBed",
			Command = f"bedtools subtract -a \"{GenomeBED}\" -b \"{CaptureBED}\" | sed -e \'s/$/\\t\\./\' > \"{NotCaptureBED}\"",
			Logger = Logger)
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.CaptureCoverage",
			Command = f"bedtools coverage -hist -sorted -g \"{RefIndex}\" -a \"{CaptureBED}\" -b \"{FinalBAM}\" | grep -P \"^all.*$\" > \"{CaptureTemp}\"",
			Logger = Logger)
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.NotCaptureCoverage",
			Command = f"bedtools coverage -hist -sorted -g \"{RefIndex}\" -a \"{NotCaptureBED}\" -b \"{FinalBAM}\" | grep -P \"^all.*$\" > \"{NotCaptureTemp}\"",
			Logger = Logger)
		# Load
		CaptureData = pandas.read_csv(CaptureTemp, sep='\t', header=None, dtype={1: int, 4: float})[[1, 4]]
		NotCaptureData = pandas.read_csv(NotCaptureTemp, sep='\t', header=None, dtype={1: int, 4: float})[[1, 4]]
	# Report Compilation
	Result = {
		"Name": [Name],
		"Capture DP>10 [%]": [CaptureData[CaptureData[1] >= 10][4].sum() * 100],
		"Capture Average": [CaptureData.apply(lambda x: x[1] * x[4], axis=1).sum()],
		"NotCapture Average": [NotCaptureData.apply(lambda x: x[1] * x[4], axis=1).sum()],
		"Enrichment Average": [None],
		"Capture DP0 [%]": [CaptureData.loc[0, 4] * 100],
		"NotCapture DP0 [%]": [NotCaptureData.loc[0, 4] * 100]
		}
	try:
		Result["Enrichment Average"] = [Result["Capture Average"][0] / Result["NotCapture Average"][0]]
	except ZeroDivisionError:
		Result["Enrichment Average"] = [float('NaN')]
	Report = "\n".join([
		f"# Scissors Coverage Stats",
		f"# Name: {Name}",
		f"# Input BAM: {FinalBAM}",
		f"# Reference: {Reference}",
		f"# Capture: {CaptureBED}",
		pandas.DataFrame(Result).to_csv(sep='\t', index=False, float_format='%.3f')])
	# Report Save
	with open(StatsTXT, 'wt') as StatsFile: StatsFile.write(Report)

def HarvestStats(
	Units: list,
	Options: dict) -> None:
	MODULE_NAME = "HarvestStats"
	
	def FormatTable(item):
		if type(item) == int: return f'{item: }'
		elif type(item) == str: return item
		elif type(item) == float: return f'{item:.4f}'.replace('.', ',')
		elif type(item) == list:
			item = [f'{i: }' if type(i) == int else i for i in item]
			item = [f'{i:.4f}'.replace('.', ',') if type(i) == float else i for i in item]
			return '; '.join(item)
	
	ReportFile = os.path.join(Options["PoolDir"], "summary.tsv")
	Data = []
	for Unit in Units:
		FileNames = GenerateFileNames(Unit, Options)
		MarkDupStats = io.StringIO(''.join(open(FileNames["DuplessMetrics"], 'rt').readlines()).split('\n\n')[1])
		MarkDupStats = pandas.read_csv(MarkDupStats, sep='\t', comment='#').apply(lambda x: x.to_list(), axis=0, result_type="reduce").to_dict()
		CoverageStats = pandas.read_csv(FileNames["CoverageStats"], sep='\t', comment='#').squeeze().to_dict()
		PrimaryStats = pandas.read_csv(FileNames["PrimaryStats"], sep='\t', names=[f"QC-passed", f"QC-failed", "index"]).set_index("index")["QC-passed"].to_dict()
		Result = {}
		LibCount = list(range(len(MarkDupStats["LIBRARY"])))
		Result["ID|"] = Unit["ID"]
		Result["Total|reads"] = int(PrimaryStats["total (QC-passed reads + QC-failed reads)"])
		Result["Mapped|% of total"] = float(PrimaryStats["mapped %"][:-1])
		Result["Secondary|% of total"] = int(PrimaryStats["secondary"]) / Result["Total|reads"] * 100
		Result["Supplementary|% of total"] = int(PrimaryStats["supplementary"]) / Result["Total|reads"] * 100
		Result["Sequenced|pairs"] = int(int(PrimaryStats["paired in sequencing"]) / 2)
		Result["Both mapped|% of sequenced"] = int(PrimaryStats["with itself and mate mapped"]) / (Result["Sequenced|pairs"] * 2) * 100
		Result["Properly paired|% of sequenced"] = float(PrimaryStats["properly paired %"][:-1])
		Result["Different chrs|% of sequenced"] = int(PrimaryStats["with mate mapped to a different chr (mapQ>=5)"]) / (Result["Sequenced|pairs"] * 2) * 100
		Result["Singletons|% of sequenced"] = float(PrimaryStats["singletons %"][:-1])
		Result["Capture DP > 10|% of capture"] = float(CoverageStats["Capture DP>10 [%]"])
		Result["Capture Average|reads by pos"] = float(CoverageStats["Capture Average"])
		Result["NotCapture Average|reads by pos"] = float(CoverageStats["NotCapture Average"])
		Result["Enrichment Average|x"] = float(CoverageStats["Enrichment Average"])
		Result["Capture DP = 0|% of capture"] = float(CoverageStats["Capture DP0 [%]"])
		Result["NotCapture DP = 0|% of non-capture"] = float(CoverageStats["NotCapture DP0 [%]"])
		Result["Libs|."] = MarkDupStats["LIBRARY"]
		Result["Total by lib|reads"] = [ int(MarkDupStats["READ_PAIRS_EXAMINED"][i]) * 2 + int(MarkDupStats["UNPAIRED_READS_EXAMINED"][i]) + int(MarkDupStats["SECONDARY_OR_SUPPLEMENTARY_RDS"][i]) + int(MarkDupStats["UNMAPPED_READS"][i]) for i in LibCount]
		Result["Mapped by lib|% of total"] = [(Result["Total by lib|reads"][i] - int(MarkDupStats["UNMAPPED_READS"][i])) / Result["Total by lib|reads"][i] * 100 for i in LibCount]
		Result["Non-primary by lib|% of total"] = [ MarkDupStats["SECONDARY_OR_SUPPLEMENTARY_RDS"][i] / Result["Total by lib|reads"][i] * 100 for i in LibCount]
		Result["Both mapped by lib|pairs"] = [int(MarkDupStats["READ_PAIRS_EXAMINED"][i]) for i in LibCount]
		Result["Singletons by lib|pairs/reads"] = [int(MarkDupStats["UNPAIRED_READS_EXAMINED"][i]) for i in LibCount]
		Result["Pair dups|% of both mapped"] = [int(MarkDupStats["READ_PAIR_DUPLICATES"][i]) / Result["Both mapped by lib|pairs"][i] * 100 for i in LibCount]
		Result["Pair optical dups|% of both mapped"] = [int(MarkDupStats["READ_PAIR_OPTICAL_DUPLICATES"][i]) / Result["Both mapped by lib|pairs"][i] * 100 for i in LibCount]
		Result["Singleton dups|% of singletons"] = [int(MarkDupStats["UNPAIRED_READ_DUPLICATES"][i]) / Result["Singletons by lib|pairs/reads"][i] * 100 for i in LibCount]
		Result["Duplication by lib|%"] = [float(x) * 100 for x in MarkDupStats["PERCENT_DUPLICATION"]]
		Result["Estimated lib size|."] = [int(MarkDupStats["ESTIMATED_LIBRARY_SIZE"][i]) for i in LibCount]
		Result = pandas.DataFrame({key: [value] for key, value in Result.items()}).transpose()
		Result = Result.applymap(FormatTable)
		Data += [ dc(Result) ]
	Data = pandas.concat(Data, axis=1).reset_index()
	Data = Data.set_index(Data["index"].apply(lambda x: x.split("|")[0]).rename(None))
	Data["index"] = Data["index"].apply(lambda x: x.split("|")[1])
	Data = Data.transpose()
	Data.to_csv(ReportFile, sep='\t', index=False)
	
def ComposeRGTag(ReadName, Sample, Library, Logger):
	
	# TODO KnV
	
	MODULE_NAME = "ComposeRGTag"
	# Regex
	LineFormats = {
		"ILLUMINA;New;Barcode": "^@(?P<Instrument>[\\w-]+):(?P<Run>\\d+):(?P<FlowCellID>[\\w-]+):(?P<Lane>\\d+):(?P<Tile>\\d+):(?P<Xpos>\\d+):(?P<Ypos>\\d+) (?P<Read>[1|2]):(?P<Filtered>[Y|N]):(?P<ControlNumber>\\d*[02468]):(?P<Barcode>[ATGCN\\+]+)$",
		"ILLUMINA;New;Sample": "^@(?P<Instrument>[\\w-]+):(?P<Run>\\d+):(?P<FlowCellID>[\\w-]+):(?P<Lane>\\d+):(?P<Tile>\\d+):(?P<Xpos>\\d+):(?P<Ypos>\\d+) (?P<Read>[1|2]):(?P<Filtered>[Y|N]):(?P<ControlNumber>\\d*[02468]):(?P<Sample>\\d+)$",
		"ILLUMINA;Old;Barcode": "^@(?P<Instrument>[\\w-]+):(?P<Lane>\\d+):(?P<Tile>\\d+):(?P<Xpos>\\d+):(?P<Ypos>\\d+)#(?P<Barcode>[ATGCN\\+]+)\\/(?P<Read>[1|2])$",
		"ILLUMINA;Old;Sample": "^@(?P<Instrument>[\\w-]+):(?P<Lane>\\d+):(?P<Tile>\\d+):(?P<Xpos>\\d+):(?P<Ypos>\\d+)#(?P<Sample>\\d+)\\/(?P<Read>[1|2])$",
		"ILLUMINA;Strange1;Barcode": "^@V(?P<Instrument>\\d+)L(?P<Lane>\\d+)C(?P<Tile>\\d+)R(?P<Coords>\\d+) (?P<Read>[1|2]):(?P<Filtered>[Y|N]):(?P<ControlNumber>\\d*[02468]):(?P<Barcode>[ATGCN\\+]+)$",
		"ILLUMINA;NCBI_1;Sample": "^@(?P<Sample>[^\\.]+).* (?P<Instrument>[\\w-]+):(?P<Run>\\d+):(?P<FlowCellID>[\\w-]+):(?P<Lane>\\d+):(?P<Tile>\\d+):(?P<Xpos>\\d+):(?P<Ypos>\\d+) length=\\d+$",
		"ILLUMINA;NCBI_2;Sample": "^@(?P<Sample>[^\\.]+).* (?P<Instrument>[\\w-]+):(?P<Lane>\\d+):(?P<Tile>\\d+):(?P<Xpos>\\d+):(?P<Ypos>\\d+) length=\\d+$",
		"ILLUMINA;NCBI_21;Sample": "^@(?P<Sample>[^\\.]+).* (?P<Instrument>[\\w-]+)\\.s_(?P<Samp>\\d+):(?P<Lane>\\d+):(?P<Tile>\\d+):(?P<Xpos>\\d+):(?P<Ypos>\\d+) length=\\d+$",
		"ILLUMINA;NCBI_ihk;Sample": "^@(?P<Sample>[^\\.]+).* SL(?P<Lane>\\d+)_R(?P<Run>\\d+)_(?P<Instrument>[\\w-]+)_.* length=\\d+$"
		}
	
	# Apply regex
	Format = [f for f in LineFormats if re.fullmatch(LineFormats[f], ReadName) is not None]
	if not Format: 
		ErrorMessage1 = f"Unknown format of line '{ReadName}'"
		Logger.error(ErrorMessage1)
		raise ValueError(ErrorMessage1)
	if len(Format) > 1:
		ErrorMessage2 = f"Ambiguous format of line '{ReadName}' (types: {', '.join(Format)})"
		Logger.error(ErrorMessage2)
		raise ValueError(ErrorMessage2)
	Fields = re.match(LineFormats[Format[0]], ReadName).groupdict()
	
	# Compose @RG
	Format = Format[0].split(";")
	PL = Format[0]
	SM = Sample
	LB = f"LIB-{Sample}-{Library}"
	if PL == "ILLUMINA":
		if Format[1] == "New": ID = f"FC{Fields['FlowCellID']}.L{Fields['Lane']}"
		if Format[1] == "Old": ID = f"FC{Fields['Instrument']}.L{Fields['Lane']}"
		if Format[1] == "Strange1": ID = f"FC{Fields['Instrument']}.L{Fields['Lane']}"
		if Format[1] == "NCBI_1": ID = f"FC{Fields['FlowCellID']}.L{Fields['Lane']}"
		if Format[1] == "NCBI_2": ID = f"FC{Fields['Instrument']}.L{Fields['Lane']}"
		if Format[1] == "NCBI_21": ID = f"FC{Fields['Instrument']}.L{Fields['Lane']}"
		if Format[1] == "NCBI_ihk": ID = f"FC{Fields['Instrument']}.L{Fields['Lane']}"
		
		PU = f"{ID}.{Fields[Format[2]]}"
	
	# Return
	return f"@RG\\tID:{ID}\\tPL:{PL}\\tPU:{PU}\\tLB:{LB}\\tSM:{SM}"

# ------======| DAEMONIC PIPELINE |======------

def DaemonicPipe(
	PipelineConfigFile: str,
	UnitsFile: str) -> None:
	MODULE_NAME = "DaemonicPipe"
	# Load data
	CurrentStage = f"{UnitsFile}.daemonic_stage"
	PipelineConfig = json.load(open(PipelineConfigFile, 'rt'))
	if os.path.exists(CurrentStage) and os.path.isfile(CurrentStage):
		Protocol = json.load(open(CurrentStage, 'rt'))
		warnings.warn(f"Resume previously interrupted pipeline from backup '{CurrentStage}'")
	else: 
		Protocol = json.load(open(UnitsFile, 'rt'))
	# Create backup
	BackupPossible = True
	try:
		SaveJSON(Protocol, CurrentStage)
	except:
		warnings.warn(f"Backup file '{CurrentStage}' cannot be created. If pipeline is interrupted, changes will not be saved")
		BackupPossible = False
	# Processing
	for Unit in Protocol["Units"]:
		# Timestamp
		StartTime = time.time()
		# Compose filenames
		FileNames = GenerateFileNames(Unit, Protocol["Options"])
		
		if Unit["Stage"] == 0:
			
			# Create dirs
			os.mkdir(FileNames["OutputDir"])
			os.mkdir(FileNames["IRs"])
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Protocol, CurrentStage)
		
		# Configure logging
		Logger = DefaultLogger(FileNames["Log"])
		
		if Unit["Stage"] == 1:
			
			# Check/add/replace @RG
			for item in Unit["Input"]:
				if item["Type"] == "fastq":
					ReadName = OpenAnyway(item["R1"], 'r', Logger).readline().decode('utf-8')[:-1]
					if item["RG"] is None: item["RG"] = ComposeRGTag(ReadName, item["Sample"], item["Library"], Logger)
				if item["Type"] == "bam": pass # TODO
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Protocol, CurrentStage)
		
		if Unit["Stage"] == 2:
			
			with tempfile.TemporaryDirectory() as TempDir:
				Shards = []
				for index, item in enumerate(Unit["Input"]):
					if item["Type"] == "fastq":
						
						TempR1, TempR2 = item["R1"], item["R2"]
						
						# Convert SOLiD
						if item["Format"] == 'solid':
							_TempR1, _TempR2 = os.path.join(TempDir, f"solid2illumina_R1_{str(index)}.fastq.gz"), None if item["R2"] is None else os.path.join(TempDir, f"solid2illumina_R2_{str(index)}.fastq.gz")
							Solid2Illumina(TempR1, _TempR1, Logger)
							if item["R2"] is not None: Solid2Illumina(TempR2, _TempR2, Logger)
							TempR1, TempR2 = _TempR1, _TempR2
						
						# Trim
						if item["Adapter"] is not None:
							ReportTXT = os.path.join(Unit['OutputDir'], f"{Unit['ID']}.InputItem{str(index)}.cutadapt.txt")
							_TempR1, _TempR2 = os.path.join(TempDir, f"adapter_R1_{str(index)}.fastq.gz"), None if item["R2"] is None else os.path.join(TempDir, f"adapter_R2_{str(index)}.fastq.gz")
							Cutadapt(TempR1, TempR2, _TempR1, _TempR2, PipelineConfig["Adapters"][item["Adapter"]], ReportTXT, PipelineConfig["Threads"], Logger)
							TempR1, TempR2 = _TempR1, _TempR2
						
						# Align
						OutputBAM = os.path.join(TempDir, f"temp_{str(index)}.bam")
						BWA(TempR1, TempR2, PipelineConfig["Reference"], item["RG"], OutputBAM, Logger, PipelineConfig["GATKCondaEnv"], Threads=PipelineConfig["Threads"])
						Shards += [ OutputBAM ]
					if item["Type"] == "bam": pass # TODO
				
				# Merge & Copy
				if len(Shards) > 1: MergeBAMs(Shards, FileNames["PrimaryBAM"], "queryname", Logger, PipelineConfig["GATKCondaEnv"])
				else: SimpleSubprocess(f"{MODULE_NAME}.CopyBAM", f"cp \"{Shards[0]}\" \"{FileNames['PrimaryBAM']}\"", Logger)
				SimpleSubprocess(f"{MODULE_NAME}.FlagStats", f"samtools flagstat -O tsv \"{FileNames['PrimaryBAM']}\" > \"{FileNames['PrimaryStats']}\"", Logger)
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Protocol, CurrentStage)
		
		if Unit["Stage"] == 3:
			
			# Mark duplicates
			MarkDuplicates(FileNames["PrimaryBAM"], FileNames["DuplessBAM"], FileNames["DuplessMetrics"], Logger, PipelineConfig["GATKCondaEnv"])
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Protocol, CurrentStage)
		
		if Unit["Stage"] == 4:
			
			# Base Recalibration
			BaseRecalibration(FileNames["DuplessBAM"], FileNames["RecalBAM"], PipelineConfig["dbSNP"], PipelineConfig["Reference"], PipelineConfig["GATK_ConfigFile"], Logger, PipelineConfig["GATKCondaEnv"], Threads=PipelineConfig["Threads"])
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Protocol, CurrentStage)
		
		if Unit["Stage"] == 5:
			
			# Coverage Stats
			CoverageStats(Unit['ID'], FileNames["RecalBAM"], FileNames["CoverageStats"], PipelineConfig["Capture"], PipelineConfig["Reference"], Logger)
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Protocol, CurrentStage)
		
		if Unit["Stage"] == 6:
			
			# Variant Calling
			HaplotypeCalling(FileNames["RecalBAM"], FileNames["VCF"], PipelineConfig["Reference"], Logger, PipelineConfig["GATKCondaEnv"], Threads=PipelineConfig["Threads"])
			
			# Timestamp
			Logger.info(f"Unit {Unit['ID']} successfully finished, summary time - %s" % (SecToTime(time.time() - StartTime)))
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Protocol, CurrentStage)
	
	HarvestStats(
		Units = Protocol["Units"],
		Options = Protocol["Options"])
	
	os.remove(CurrentStage)
