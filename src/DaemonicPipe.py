from src.SharedFunctions import *


## ------======| PIPELINE STAGES |======------

def Solid2Illumina(InputFQ, OutputFQ, Logger):
	
	MODULE_NAME = "Solid2Illumina"
	
	# Logging
	Logger.info(f"Input FASTQ: {InputFQ}")
	Logger.info(f"Output FASTQ: {OutputFQ}")
	
	# Processing
	SimpleSubprocess(f"{MODULE_NAME}.Convert", f"cutadapt -c --format=sra-fastq --bwa --action=none -o \"{OutputFQ}\" \"{InputFQ}\"", Logger)
	
def Cutadapt(InputR1, InputR2, OutputR1, OutputR2, Adapter, ReportTXT, Threads, Logger):
	
	MODULE_NAME = "Cutadapt"
	
	# Logging
	if InputR2 is not None:
		Logger.info(f"Mode: Paired-end")
		Logger.info(f"Input FASTQ [R1]: {InputR1}")
		Logger.info(f"Input FASTQ [R2]: {InputR2}")
		Logger.info(f"Output FASTQ [R1]: {OutputR1}")
		Logger.info(f"Output FASTQ [R2]: {OutputR2}")
	else:
		Logger.info(f"Mode: Single-end")
		Logger.info(f"Input FASTQ: {InputR1}")
		Logger.info(f"Output FASTQ: {OutputR1}")
	
	Logger.info(f"Adapter: {Adapter['Name']}")
	
	# Processing
	if InputR2 is not None: SimpleSubprocess(f"{MODULE_NAME}.Trim", f"cutadapt -j {str(Threads)} -e 0.2 -m 8 -a {Adapter['R1']} -A {Adapter['R2']} -o \"{OutputR1}\" -p \"{OutputR2}\" \"{InputR1}\" \"{InputR2}\" > \"{ReportTXT}\"", Logger)
	else: SimpleSubprocess(f"{MODULE_NAME}.Trim", f"cutadapt -j {str(Threads)} -e 0.2 -m 8 -a {Adapter['R1']} -o \"{OutputR1}\" \"{InputR1}\" > \"{ReportTXT}\"", Logger)

def BamMetrics(BamMetricsFile):
	with open(f"/dev/datasets/FairWind/_results/The_New_96/MISEQ_2020-{item}/MISEQ_2020-{item}.primary_stats.txt", 'rt') as f: buff = '::'.join([item[:-1] for item in f])
	primary_stats = re.match("^(?P<TotalPassed>\\d+) \\+ (?P<TotalFailed>\\d+) in total \\(QC-passed reads \\+ QC-failed reads\\)::(?P<SecondaryPassed>\\d+) \\+ (?P<SecondaryFailed>\\d+) secondary::(?P<SupplementaryPassed>\\d+) \\+ (?P<SupplementaryFailed>\\d+) supplementary::(?P<DuplicatesPassed>\\d+) \\+ (?P<DuplicatesFailed>\\d+) duplicates::(?P<MappedPassed>\\d+) \\+ (?P<MappedFailed>\\d+) mapped.*::(?P<PairedPassed>\\d+) \\+ (?P<PairedFailed>\\d+) paired in sequencing::(?P<Read1Passed>\\d+) \\+ (?P<Read1Failed>\\d+) read1::(?P<Read2Passed>\\d+) \\+ (?P<Read2Failed>\\d+) read2::(?P<ProperlyPairedPassed>\\d+) \\+ (?P<ProperlyPairedFailed>\\d+) properly paired.*::(?P<BothMappedPassed>\\d+) \\+ (?P<BothMappedFailed>\\d+) with itself and mate mapped::(?P<SingletonsPassed>\\d+) \\+ (?P<SingletonsFailed>\\d+) singletons.*::(?P<MateOnDifferentChrPassed>\\d+) \\+ (?P<MateOnDifferentChrFailed>\\d+) with mate mapped to a different chr::(?P<MateOnDifferentChrMapQ5Passed>\\d+) \\+ (?P<MateOnDifferentChrMapQ5Failed>\\d+) with mate mapped to a different chr.*$", buff).groupdict()
	pass # TODO

def PercentMarkDupMetrics(MDMetricsFile):
	for item in pandas.read_csv(MDMetricsFile, sep='\t', comment='#', chunksize=1): return item["PERCENT_DUPLICATION"][0]

def CoverageStats(Name, FinalBAM, StatsTXT, CaptureBED, NotCaptureBED, Reference, Logger):
	
	MODULE_NAME = "CoverageStats"
	
	# Logging
	Logger.info(f"BAM File: {FinalBAM}")
	Logger.info(f"Reference: {Reference}")
	Logger.info(f"Capture BED: {CaptureBED}")
	Logger.info(f"NOT Capture BED: {NotCaptureBED}")
	
	# Options
	RefIndex = f"{Reference}.fai"
	MedianFunc = lambda List: [num for num, item in enumerate(List) if sum(List[num:]) >= 0.5][-1]
	AverageFunc = lambda List: sum([num * item for num, item in enumerate(List)])
	
	with tempfile.TemporaryDirectory() as TempDir:
		
		CaptureTemp = os.path.join(TempDir, "capture.csv")
		NotCaptureTemp = os.path.join(TempDir, "not_capture.csv")
		
		# Coverage
		SimpleSubprocess(f"{MODULE_NAME}.CaptureCoverage", f"bedtools coverage -hist -sorted -g \"{RefIndex}\" -a \"{CaptureBED}\" -b \"{FinalBAM}\" | grep -P \"^all.*$\" > \"{CaptureTemp}\"", Logger)
		SimpleSubprocess(f"{MODULE_NAME}.NotCaptureCoverage", f"bedtools coverage -hist -sorted -g \"{RefIndex}\" -a \"{NotCaptureBED}\" -b \"{FinalBAM}\" | grep -P \"^all.*$\" > \"{NotCaptureTemp}\"", Logger)
		
		# Load
		CaptureData = pandas.read_csv(CaptureTemp, sep='\t', header=None, dtype={4: float})[4].to_list()
		NotCaptureData = pandas.read_csv(NotCaptureTemp, sep='\t', header=None, dtype={4: float})[4].to_list()
		
	# Report Compilation
	ResultIndex = ["Name", "Capture DP>10 [%]", "Capture Average", "NotCapture Average", "Capture Median", "NotCapture Median", "Enrichment Average", "Enrichment Median", "Capture DP0 [%]", "NotCapture DP0 [%]"]
	Result = pandas.Series([None] * len(ResultIndex), index=ResultIndex)
	
	Result["Name"] = Name
	Result["Capture DP>10 [%]"] = sum(CaptureData[10:]) * 100
	Result["Capture DP0 [%]"] = CaptureData[0] * 100
	Result["NotCapture DP0 [%]"] = NotCaptureData[0] * 100
	Result["Capture Median"] = MedianFunc(CaptureData)
	Result["NotCapture Median"] = MedianFunc(NotCaptureData)
	Result["Enrichment Median"] = Result["Capture Median"] / (Result["NotCapture Median"] + 1)
	Result["Capture Average"] = AverageFunc(CaptureData)
	Result["NotCapture Average"] = AverageFunc(NotCaptureData)
	try:
		Result["Enrichment Average"] = Result["Capture Average"] / Result["NotCapture Average"]
	except ZeroDivisionError:
		Result["Enrichment Average"] = float('NaN')
	
	Report = "# Scissors Coverage Stats\n# Name: " + Name + "\n# Input BAM: " + FinalBAM + "\n# Reference: " + Reference + "\n# Capture: " + CaptureBED + "\n# NOT Capture: " + NotCaptureBED + "\n" + Result.to_frame().transpose().to_csv(sep='\t', index=False, float_format='%.2f')
	with open(StatsTXT, 'wt') as StatsFile: StatsFile.write(Report)

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

def FastQC(InputFastQ, OutputHTML, Logger, Size=0, Threads=cpu_count()):
	
	MODULE_NAME = "FastQC"
	
	# Logging
	Logger.info(f"FastQ: {InputFastQ}")
	Logger.info(f"Report: {OutputHTML}")
	
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		AnalyzeFilename = InputFastQ
		if Size != 0:
			Logger.info(f"Subsample: {str(Size)}")
			SampleFilename = os.path.join(TempDir, "sample.fastq.gz")
			SimpleSubprocess(f"{MODULE_NAME}.Sampling", f"zcat -q \"{InputFastQ}\" | head -{str(Size * 4)} | gzip -c > \"{SampleFilename}\"", Logger)
			AnalyzeFilename = SampleFilename
		SimpleSubprocess(f"{MODULE_NAME}.Analysis", f"fastqc -o \"{TempDir}\" -t {str(Threads)} \"{AnalyzeFilename}\"", Logger)
		HTMLTemp = glob(os.path.join(TempDir, "*.html"))
		if len(HTMLTemp) != 1:
			ErrorMessage = f"Error processing file '{InputFastQ}'"
			Logger.error(ErrorMessage)
			raise RuntimeError(ErrorMessage)
		SimpleSubprocess(f"{MODULE_NAME}.Move", f"cp \"{HTMLTemp[0]}\" \"{OutputHTML}\"", Logger)

def BWA(InputR1, InputR2, Reference, RGHeader, OutputBAM, Logger, Env, Threads=cpu_count()):
	
	MODULE_NAME = "BWA"
	
	# Logging
	if InputR2 is None:
		Logger.info(f"Reads: {InputR1}")
	else:
		Logger.info(f"R1: {InputR1}")
		Logger.info(f"R2: {InputR2}")
	Logger.info(f"Output BAM: {OutputBAM}")
	Logger.info(f"Reference: {Reference}")
	Logger.info(f"RG Header: {RGHeader}")
	
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		InputR2 = "" if (InputR2 is None) else f"\"{InputR2}\""
		SimpleSubprocess(f"{MODULE_NAME}.AlignAndSort", f"bwa mem -R \"{RGHeader}\" -t {str(Threads)} -v 1 \"{Reference}\" \"{InputR1}\" {InputR2} | gatk SortSam --VERBOSITY ERROR --TMP_DIR \"{TempDir}\" -SO queryname -I \"/dev/stdin\" -O \"{OutputBAM}\"", Logger, CheckPipefail=True, Env=Env)

def MergeBAMs(BAMs, OutputBAM, SortOrder, Logger, Env):
	
	MODULE_NAME = "MergeBAMs"
	
	# Processing
	SimpleSubprocess(f"{MODULE_NAME}.Merge", f"gatk MergeSamFiles --USE_THREADING true -SO {SortOrder} {MultipleTags('-I', BAMs)} -O \"{OutputBAM}\"", Logger, Env=Env)

def MarkDuplicates(InputBAM, OutputBAM, MetricsTXT, Logger, Env):
	
	MODULE_NAME = "MarkDuplicates"
	
	# Logging
	Logger.info(f"Input: {InputBAM}")
	Logger.info(f"Output: {OutputBAM}")
	Logger.info(f"Metrics: {MetricsTXT}")
	
	# Options
	JavaOptions = "-XX:+UseParallelGC -XX:ParallelGCThreads=2" # -Xmx12G removed, can't see the sense of memory limits if the process is single-threaded.
	
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		SimpleSubprocess(f"{MODULE_NAME}.RemoveAndSort", f"gatk --java-options \"{JavaOptions}\" MarkDuplicates --REMOVE_DUPLICATES true --VERBOSITY ERROR --ASSUME_SORT_ORDER queryname --TMP_DIR \"{TempDir}\" -M \"{MetricsTXT}\" -I \"{InputBAM}\" -O \"/dev/stdout\" | gatk SortSam --VERBOSITY ERROR --TMP_DIR \"{TempDir}\" -SO coordinate -I \"/dev/stdin\" -O \"{OutputBAM}\"", Logger, CheckPipefail=True, Env=Env)

def ContigBaseRecalibration(Contig, InputBAM, TempDir, dbSNP, Reference, GATK_ConfigFile, Logger, Env):
	
	MODULE_NAME = "ContigBaseRecalibration"
	
	# Options
	FiltersComparison = ["MappedReadFilter", "MappingQualityAvailableReadFilter", "MappingQualityNotZeroReadFilter", "NotDuplicateReadFilter", "NotSecondaryAlignmentReadFilter", "PassesVendorQualityCheckReadFilter"] 
	# Fucking bug. BaseRecalibrator and ApplyBQSR filter reads differently, so ApplyBQSR fucks up every time it can't find read group.
	BaseRecalibrator_JavaOptions = "-Xmx3G -XX:+UseParallelGC -XX:ParallelGCThreads=2"
	ApplyBQSR_JavaOptions = "-Xmx3G"
	BQSRTable = os.path.join(TempDir, f"bqsr_table_{Contig}.tsv")
	OutputBAM = os.path.join(TempDir, f"output_{Contig}.bam")
	
	# Processing
	SimpleSubprocess(f"{MODULE_NAME}.MakeTable[{Contig}]", f"gatk --java-options \"{BaseRecalibrator_JavaOptions}\" BaseRecalibrator --gatk-config-file \"{GATK_ConfigFile}\" --tmp-dir \"{TempDir}\" -L {Contig} -I \"{InputBAM}\" --known-sites \"{dbSNP}\" -O \"{BQSRTable}\" -R \"{Reference}\"", Logger, Env=Env)
	SimpleSubprocess(f"{MODULE_NAME}.Apply[{Contig}]", f"gatk --java-options \"{ApplyBQSR_JavaOptions}\" ApplyBQSR --gatk-config-file \"{GATK_ConfigFile}\" {MultipleTags('-RF', FiltersComparison, Quoted=False)} --tmp-dir \"{TempDir}\" -OBI false -L {Contig} -bqsr \"{BQSRTable}\" -I \"{InputBAM}\" -O \"{OutputBAM}\"", Logger, Env=Env)
	
	# Return
	return OutputBAM

def BaseRecalibration(InputBAM, OutputBAM, dbSNP, Reference, GATK_ConfigFile, Logger, Env, Threads=cpu_count()):
	
	MODULE_NAME = "BaseRecalibration"
	
	# Logging
	Logger.info(f"Input: {InputBAM}")
	Logger.info(f"Output: {OutputBAM}")
	Logger.info(f"Known sites: {dbSNP}")
	Logger.info(f"Reference: {Reference}")
	
	# Options
	Contigs = [item['SN'] for item in GetContigs(InputBAM)]
	
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		SimpleSubprocess(f"{MODULE_NAME}.PreIndex", f"gatk BuildBamIndex -I \"{InputBAM}\"", Logger, Env=Env)
		with Threading("ContigBaseRecalibration", Logger, Threads) as pool: Shards = pool.map(functools.partial(ContigBaseRecalibration, InputBAM=InputBAM, TempDir=TempDir, dbSNP=dbSNP, Reference=Reference, GATK_ConfigFile=GATK_ConfigFile, Logger=Logger, Env=Env), Contigs)
		SimpleSubprocess(f"{MODULE_NAME}.Merge", f"gatk MergeSamFiles --USE_THREADING true -SO coordinate {MultipleTags('-I', Shards)} -O \"{OutputBAM}\"", Logger, Env=Env)
		SimpleSubprocess(f"{MODULE_NAME}.PostIndex", f"gatk BuildBamIndex -I \"{OutputBAM}\"", Logger, Env=Env)

def ContigHaplotypeCalling(Contig, InputBAM, TempDir, Reference, Logger, Env):
	
	MODULE_NAME = "ContigHaplotypeCalling"
	
	# Options
	JavaOptions = "-Xmx3G"
	OutputVCF = os.path.join(TempDir, f"output_{Contig}.vcf")
	
	# Processing
	SimpleSubprocess(f"{MODULE_NAME}.Calling[{Contig}]", f"gatk --java-options \"{JavaOptions}\" HaplotypeCaller --native-pair-hmm-threads 2 -OVI false --dont-use-soft-clipped-bases true -L {Contig} -I \"{InputBAM}\" -O \"{OutputVCF}\" -R \"{Reference}\"", Logger, Env=Env)
	
	# Return
	return OutputVCF

def HaplotypeCalling(InputBAM, OutputVCF, Reference, Logger, Env, Threads=cpu_count()):
	
	MODULE_NAME = "HaplotypeCalling"
	
	# Logging
	Logger.info(f"Input BAM: {InputBAM}")
	Logger.info(f"Output VCF: {OutputVCF}")
	Logger.info(f"Reference: {Reference}")
	
	# Options
	Contigs = [item['SN'] for item in GetContigs(InputBAM)]
	
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		with Threading("ContigHaplotypeCalling", Logger, Threads) as pool: Shards = pool.map(functools.partial(ContigHaplotypeCalling, InputBAM=InputBAM, TempDir=TempDir, Reference=Reference, Logger=Logger, Env=Env), Contigs)
		SimpleSubprocess(f"{MODULE_NAME}.Merge", f"gatk MergeVcfs {MultipleTags('-I', Shards)} -O \"{OutputVCF}\"", Logger, Env=Env)

def VariantRecalibration(InputBAM, InputVCF, AnnotatedVCF, OutputVCF, Reference, Resources, Tranches, Logger, Env, DisableAVX=False):
	
	MODULE_NAME = "VariantRecalibration"
	
	# Logging
	Logger.info(f"Input BAM: {InputBAM}")
	Logger.info(f"Input VCF: {InputVCF}")
	Logger.info(f"Output VCF: {OutputVCF}")
	Logger.info(f"Reference: {Reference}")
	for index, item in enumerate(Resources): Logger.info(f"Resource #{str(index + 1)}: {item}")
	Logger.info(f"SNP Tranches: {'; '.join([str(item) for item in Tranches['SNP']])}")
	Logger.info(f"INDEL Tranches: {'; '.join([str(item) for item in Tranches['INDEL']])}")
	
	# Processing
	SimpleSubprocess(f"{MODULE_NAME}.ScoreVariants", f"gatk CNNScoreVariants --disable-avx-check {'true' if DisableAVX else 'false'} -I \"{InputBAM}\" -V \"{InputVCF}\" -R \"{Reference}\" -O \"{AnnotatedVCF}\" -tensor-type read_tensor", Logger, Env=Env)
	SimpleSubprocess(f"{MODULE_NAME}.FilterVariants", f"gatk FilterVariantTranches -V \"{AnnotatedVCF}\" {MultipleTags('--resource', Resources)} --info-key CNN_2D {MultipleTags('--snp-tranche', Tranches['SNP'], Quoted=False)} {MultipleTags('--indel-tranche', Tranches['INDEL'], Quoted=False)} -O \"{OutputVCF}\"", Logger, Env=Env)


# ------======| DAEMONIC PIPELINE |======------

def DaemonicPipe(PipelineConfigFile, UnitsFile):
	
	MODULE_NAME = "DaemonicPipe"
	
	# Load data
	CurrentStage = f"{UnitsFile}.daemonic_stage"
	PipelineConfig = json.load(open(PipelineConfigFile, 'rt'))
	if os.path.exists(CurrentStage) and os.path.isfile(CurrentStage): 
		Units = json.load(open(CurrentStage, 'rt'))
		warnings.warn(f"Resume previously interrupted pipeline from backup '{CurrentStage}'")
	else: Units = json.load(open(UnitsFile, 'rt'))
	
	# Create backup
	BackupPossible = True
	try:
		SaveJSON(Units, CurrentStage)
	except:
		warnings.warn(f"Backup file '{CurrentStage}' cannot be created. If pipeline is interrupted, changes will not be saved")
		BackupPossible = False
	
	# Processing
	for Unit in Units:
		
		# Timestamp
		StartTime = time.time()
		
		# Compose filenames
		IRs = os.path.join(Unit['OutputDir'], "IRs")
		FileNames = {
			"Log": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.pipeline_log.txt"),
			"PrimaryBAM": os.path.join(IRs, f"{Unit['ID']}.primary.bam"),
			"PrimaryStats": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.primary_stats.txt"),
			"DuplessBAM": os.path.join(IRs, f"{Unit['ID']}.dupless.bam"),
			"DuplessMetrics": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.md_metrics.txt"),
			"RecalBAM": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.final.bam"),
			"VCF": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.unfiltered.vcf"),
			"CoverageStats": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.coverage.txt")
			}
		
		if Unit["Stage"] == 0:
			
			# Create dirs
			os.mkdir(Unit['OutputDir'])
			os.mkdir(IRs)
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Units, CurrentStage)
		
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
			if BackupPossible: SaveJSON(Units, CurrentStage)
		
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
				SimpleSubprocess(f"{MODULE_NAME}.FlagStats", f"samtools flagstat \"{FileNames['PrimaryBAM']}\" > \"{FileNames['PrimaryStats']}\"", Logger)
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Units, CurrentStage)
		
		if Unit["Stage"] == 3:
			
			# Mark duplicates
			MarkDuplicates(FileNames["PrimaryBAM"], FileNames["DuplessBAM"], FileNames["DuplessMetrics"], Logger, PipelineConfig["GATKCondaEnv"])
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Units, CurrentStage)
		
		if Unit["Stage"] == 4:
			
			# Base Recalibration
			BaseRecalibration(FileNames["DuplessBAM"], FileNames["RecalBAM"], PipelineConfig["dbSNP"], PipelineConfig["Reference"], PipelineConfig["GATK_ConfigFile"], Logger, PipelineConfig["GATKCondaEnv"], Threads=PipelineConfig["Threads"])
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Units, CurrentStage)
		
		if Unit["Stage"] == 5:
			
			# Coverage Stats
			CoverageStats(Unit['ID'], FileNames["RecalBAM"], FileNames["CoverageStats"], PipelineConfig["Capture"], PipelineConfig["NotCapture"], PipelineConfig["Reference"], Logger)
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Units, CurrentStage)
		
		if Unit["Stage"] == 6:
			
			# Variant Calling
			HaplotypeCalling(FileNames["RecalBAM"], FileNames["VCF"], PipelineConfig["Reference"], Logger, PipelineConfig["GATKCondaEnv"], Threads=PipelineConfig["Threads"])
			
			# Timestamp
			Logger.info(f"Unit {Unit['ID']} successfully finished, summary time - %s" % (SecToTime(time.time() - StartTime)))
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Units, CurrentStage)
	
	os.remove(CurrentStage)
