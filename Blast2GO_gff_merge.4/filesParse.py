import pandas as pd
from BCBio import GFF
from goCounts import outputReview


@outputReview
def b2GoTable(file, gene_col_index, GO_term_col_index):
	"""
	Method to omicxbox_table.txt eggnog_go_slim.txt and eggnog.txt format
	Args:
	     file: file path
		 gene_col_index: input gene column index
		 GO_term_col_index: input GO terms column index
	Returns: 
	     Dictionary of Key: gene, value: list of GO terms(GO:0016977)
		 3 files format for Not matching retuns ['']
	"""
	output_dic= {} 
	table = pd.read_csv(file, sep='\t', keep_default_na=False)
	for i in range(len(table)):
		go_terms = table.iloc[i, GO_term_col_index].split('; ')
		go_terms = [i[2:] for i in go_terms]
		if "" in go_terms:
			go_terms.remove("")
		output_dic[table.iloc[i, gene_col_index]] = go_terms
	return output_dic


@outputReview
def b2GoAnnot(file):
	"""
	Method to blast2go.annot format
	Args:
	     file: file path
	Returns: 
	     Dictionary of Key: gene, value: list of GO terms(GO:0016977)
		 All inputs have a value, should be not gene without matching go term
	"""
	output_dic= {} 
	file = open(file)
	lines = file.readlines()
	for row in lines:
		row = row.replace("\n", "")
		row = row.rstrip("\t")
		row = row.split("\t")
		gene = row[0]
		go_term = row[1]
		if gene in output_dic.keys():
			output_dic[gene]=output_dic[gene] + [go_term]
		else:		
			output_dic[gene]=[go_term]

	return output_dic

@outputReview
def b2GoGffExport(file):
	"""
	Method to gff_export.gff format
	Args:
	     file: file path
	Returns: 
	     Dictionary of Key: gene, value: list of GO terms(GO:0016977)
		 Not matching retuns ['']
	"""
	output_dic = {}
	for scaffold in GFF.parse(file):
		for gene in scaffold.features:
			if gene.type == "CDS":
				qualifier = gene.qualifiers
				id_gff = qualifier.get("ID", None)[0]
				GO_gff = qualifier.get("Ontology_id", [""])
				if "true" in GO_gff:
					GO_gff=[""]
				if "" in GO_gff:
					GO_gff.remove("")
				output_dic[id_gff] = GO_gff
	return output_dic

@outputReview
def b2GOIprcscn(file):
	"""
	Method to iprscan_Tatro_V3_annotations.txt format
	Args:
	     file: file path
	Returns: 
	     Dictionary of Key: gene, value: list of GO terms(GO:0016977)
		 Not matching retuns ['-']
	"""
	output_dic = {}
	table=pd.read_csv(file,sep='\t',keep_default_na=False,index_col=False)
	for i in range(len(table)):
		go_terms=table.iloc[i, 9].split(",")
		if "-" in go_terms:
			go_terms.remove("-")
		output_dic[table.iloc[i,0]]=go_terms	
	return output_dic

