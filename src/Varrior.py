from SharedFunctions import *
import vcf

## ------======| TABIX OBJECT CLASS |======------

class LovelyTabix:
	
	def __init__(self,
			Name: str,
			TabixFile: str,
			Type: str,
			Scheme: str,
			Logger: logging.Logger):
		Types = ['coord', 'region']
		if Type not in Types: raise ValueError(f"Unknown db type: '{Type}'. Valid values: {str(Types)}")
		self.__Name, self.__TabixFile, self.__Type, self.__Scheme, self.__Logger = Name, TabixFile, Type, Scheme, Logger
	
	@classmethod
	def ParseByScheme(Line, Scheme):
		if Scheme == 'vcf':
			vcf.model._Record('\t'.join(Line))
	
	@classmethod
	def NullVcfReader(Class): return vcf.Reader(io.StringIO('.'))
	
	@classmethod
	def ParseVcfLine(Class, Line, VcfParser):
		ParsedLine = {}
		ParsedLine["CHROM"] = Line[0]
		ParsedLine["POS"] = int(Line[1])
		ParsedLine["ID"] = Line[2] if Line[2] != '.' else None
		ParsedLine["REF"] = Line[3]
		ParsedLine["ALT"] = list(map(VcfParser._parse_alt, Line[4].split(',')))
		ParsedLine["INFO"] = VcfParser._parse_info(Line[7])
		return ParsedLine
	
	@classmethod
	def SeekerThread(Class,
			Chunk: list,
			Name: str,
			TabixFile: str,
			Scheme: str,
			Logger: logging.Logger) -> list:
		MODULE_NAME = "SeekerThread"
		ThreadNum, Chunk = Chunk
		Chunk = [(str(item[0]), int(item[1]), int(item[2])) for item in Chunk]
		Seeker = tabix.open(TabixFile)
		VcfParser = LovelyTabix.NullVcfReader()
		Result = {}
		for Query in Chunk:
			Result[Query] = list(Seeker.query(Query[0], Query[1], Query[2]))
			if Scheme == 'vcf': Result[Query] = [LovelyTabix.ParseVcfLine(Line, VcfParser) for Line in Result[Query]]
		return Result
	
	def Seeker(self,
			PositionsList: list,
			Threads: int = cpu_count()) -> dict:
		MODULE_NAME = "Seeker"
		Chunks = [Chunk.tolist() for Chunk in numpy.array_split(PositionsList, Threads)]
		with Threading(f"{self.__Name}.{MODULE_NAME}", self.__Logger, Threads) as ThePool:
			Chunks = ThePool.map(functools.partial(LovelyTabix.SeekerThread, Name=self.__Name, TabixFile=self.__TabixFile, Scheme=self.__Scheme, Logger=self.__Logger), enumerate(Chunks))
			Result = {}
			for Index in range(Threads): Result = {**Result, **Chunks[Index]}
		return Result

## ------======| MAIN |======------

if __name__ == "__main__":
	Logger = DefaultLogger("/dev/null")
	to = LovelyTabix("dbSNP", "/dev/datasets/FairWind/_db/dbsnp/00-All.vcf.gz", 'coord', 'vcf', Logger)
	data = pandas.read_csv("/dev/datasets/FairWind/_results/DOLG/dcDOL4/dcDOL4.unfiltered.vcf", sep='\t', comment='#', header=None, dtype=str, nrows=1000)
	variants = data[[0,1]].apply(lambda x: (x[0][3:], int(x[1]) - 1, int(x[1])), axis=1).to_list()
	print({index: item for index, item in to.Seeker(variants).items() if item})
	
