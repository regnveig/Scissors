import quicksect
import pandas

class Regenome:
	
	def __init__(self, dbname=None, filename=None, dbtype="default"):
		self.__tree = dict()
		self.__dbtypes = ["default", "ucsc"]
		if dbtype == self.__dbtypes[0]: pass
		elif dbtype == self.__dbtypes[1]:
			assert (type(dbname) == str) and (type(filename) == str), f"dbname and filename must be string type"
			chrom, start, end = "#chrom", "chromStart", "chromEnd"
			data = pandas.read_csv(filename, sep="\t")
			for index, line in data.iterrows():
				annotation = dict(line[[item for item in line.index.to_list() if item not in [chrom, start, end]]])
				annotation = {f"{dbname}.{str(key)}": val for key, val in annotation.items()}
				if not annotation: annotation = {dbname: True}
				self.add(line[chrom], line[start], line[end], annotation=annotation)
		else: raise ValueError(f"Incorrect dbname (available names: {', '.join(self.__dbtypes)})")
	
	def __base_assertion(self, chrom, start, end):
		assert type(chrom) == str, 			f"Chrom must be string type"
		assert type(start) == int, 			f"Start position must be integer type"
		assert type(end) == int, 			f"End position must be integer type"
		assert start >= 0, 					f"Start position must be positive or zero"
		assert end >= 0, 					f"End position must be positive or zero"
		assert start <= end, 				f"Start position must be less than or equal to end position"
	
	def add(self, chrom, start, end, annotation):
		self.__base_assertion(chrom, start, end)
		assert type(annotation) == dict, 	f"Annotation must be dict type"
		try:
			self.__tree[chrom]
		except KeyError:
			self.__tree[chrom] = quicksect.IntervalTree()
		self.__tree[chrom].insert(quicksect.Interval(start, end, data=annotation))
	
	def annotation(self, chrom, start, end):
		self.__base_assertion(chrom, start, end)
		try:
			self.__tree[chrom]
		except KeyError:
			return list()
		return [item.data for item in self.__tree[chrom].search(start, end)]

class Regenome2:
	
	def __init__(self, dbname, filename, dbtype):
		self.__tree = dict()
		self.__dbtypes = ["ucsc"]
		if dbtype == self.__dbtypes[0]:
			assert (type(dbname) == str) and (type(filename) == str), f"dbname and filename must be string type"
			chrom, start, end = "#chrom", "chromStart", "chromEnd"
			data = pandas.read_csv(filename, sep="\t")
			data.index = data.index.map(lambda x: "item_" + str(x))
			for c in list(set(data["#chrom"].to_list())):
				t = data[data["#chrom"] == c][["chromStart", "chromEnd"]]
				t["chromStart"] = t["chromStart"].apply(lambda x: x - 0.1)
				t["chromEnd"] = t["chromEnd"].apply(lambda x: x + 0.1)
				t = pandas.DataFrame(t.stack()).reset_index().rename(columns={"level_0": "item", "level_1": "type", 0: "pos"}).sort_values(by="pos")
				t["index"] = t.apply(lambda x: [data.loc[key,:].to_dict() for key, value in t[t["pos"] <= x["pos"]]["item"].value_counts().to_dict().items() if value == 1], axis=1)
				print(t)
				self.__tree[c] = t
		else: raise ValueError(f"Incorrect dbname (available names: {', '.join(self.__dbtypes)})")
	
	def __base_assertion(self, chrom, start, end):
		assert type(chrom) == str, 			f"Chrom must be string type"
		assert type(start) == int, 			f"Start position must be integer type"
		assert type(end) == int, 			f"End position must be integer type"
		assert start >= 0, 					f"Start position must be positive or zero"
		assert end >= 0, 					f"End position must be positive or zero"
		assert start <= end, 				f"Start position must be less than or equal to end position"
	
	def add(self, chrom, start, end, annotation):
		self.__base_assertion(chrom, start, end)
		assert type(annotation) == dict, 	f"Annotation must be dict type"
		try:
			self.__tree[chrom]
		except KeyError:
			self.__tree[chrom] = quicksect.IntervalTree()
		self.__tree[chrom].insert(quicksect.Interval(start, end, data=annotation))
	
	def annotation(self, chrom, start, end):
		self.__base_assertion(chrom, start, end)
		try:
			self.__tree[chrom]
		except KeyError:
			return list()
		return [item.data for item in self.__tree[chrom].search(start, end)]
	
Regenome2("u", "/dev/datasets/FairWind/.cloud/core/labjournal/tools/Scissors/db/ProblematicRegions/GIAB_GenotypeConflict.canonical.tsv", "ucsc")
