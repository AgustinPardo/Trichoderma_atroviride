import pandas as pd
from BCBio import GFF

def b2GoAnnot(file):
	"""
	Method to blast2go.annot format
	Args:
	     file: file path
	Returns: 
	     {gene: {Ontology_term: value}}
		 All inputs have a value, should be not gene without matching go term
	"""
	output_dic = {}
	file = open(file)
	lines = file.readlines()
	for row in lines:
		row = row.replace("\n", "")
		row = row.rstrip("\t")
		row = row.split("\t")
		gene = row[0]
		ontology_term = row[1]
		if gene in output_dic.keys():
			output_dic[gene]["Ontology_term"] = output_dic[gene]["Ontology_term"] + [ontology_term]
		else:
			output_dic[gene] = {"Ontology_term": ontology_term}

	return output_dic


def eggnogTxt(file):
	"""
	Method to eggnog.txt format
	Args:
	     file: file path
	Returns: 
	    {gene: {product: value, note: value, Ontology_term: [value]}}
		Not matching retuns ''
	"""
	output_dic = {}
	table = pd.read_csv(file, sep='\t', keep_default_na=False)

	for i in range(len(table)):
		product_term = table.iloc[i, 2]		
		note_term = table.iloc[i, 7]
		ontology_term = table.iloc[i, 9].split('; ')
		ontology_term = [i[2:] for i in ontology_term]

		if "" in ontology_term:
			ontology_term.remove("")

		output_dic[table.iloc[i, 1]] = {"product": product_term, "note": note_term, "Ontology_term": ontology_term}
	return output_dic



def omicsboxTableTxt(file):
	"""
	Method to omicsbox_table.txt format
	Args:
	     file: file path
	Returns: 
	    {gene: {"Dbxref": {"InterPro": value}, "Ontology_term": value))}}
		Not matching retuns ''
	"""
	output_dic = {}
	table = pd.read_csv(file, sep='\t', keep_default_na=False)

	for i in range(len(table)):
		dbxref_term = table.iloc[i, 13]	

		dbxref_term = dbxref_term.split(" ")
		interPro_term = [term for term in dbxref_term if "IPR" in term]

		ontology_term = table.iloc[i, 9].split('; ')
		ontology_term = [i[2:] for i in ontology_term]

		ontology_term.remove("") if "" in ontology_term else None		

		interPro_go_ids_term=table.iloc[i, 14].split('; ')
		interPro_go_ids_term = [i[2:] for i in interPro_go_ids_term]

		interPro_go_ids_term.remove("no GO terms") if "no GO terms" in interPro_go_ids_term else None
		interPro_go_ids_term.remove(" GO terms") if " GO terms" in interPro_go_ids_term else None
		interPro_go_ids_term.remove(" IPS match") if " IPS match" in interPro_go_ids_term else None
		
		output_dic[table.iloc[i, 2]] = {"Dbxref": {"InterPro": interPro_term}, "Ontology_term": list(set(ontology_term + interPro_go_ids_term))}
	return output_dic



def omicsboxTableDiamondBlastTxt(file):
	"""
	Method to omicsbox_table_DiamondBlast.txt format
	Args:
	     file: file path
	Returns:
			{gene: {"Dbxref": {"InterPro": value}, "note": {"EC": value}, "Ontology_term": value}}
 
		Not matching retuns ''
	"""
	output_dic = {}
	table = pd.read_csv(file, sep='\t', keep_default_na=False)

	for i in range(len(table)):
		dbxref_term = table.iloc[i, 13]	

		dbxref_term = dbxref_term.split(" ")
		interPro_term = [term for term in dbxref_term if "IPR" in term]

		enzyme_code_term = table.iloc[i, 11].split('; ')
		enzyme_code_term = [i for i in enzyme_code_term]

		ontology_term = table.iloc[i, 9].split('; ')
		ontology_term = [i[2:] for i in ontology_term]

		ontology_term.remove("") if "" in ontology_term else None		

		interPro_go_ids_term=table.iloc[i, 14].split('; ')
		interPro_go_ids_term = [i[2:] for i in interPro_go_ids_term]

		interPro_go_ids_term.remove("no GO terms") if "no GO terms" in interPro_go_ids_term else None
		interPro_go_ids_term.remove(" GO terms") if " GO terms" in interPro_go_ids_term else None
		interPro_go_ids_term.remove(" IPS match") if " IPS match" in interPro_go_ids_term else None
		
		output_dic[table.iloc[i, 2]] = {"Dbxref": {"InterPro": interPro_term}, "note": {"EC": enzyme_code_term}, "Ontology_term": list(set(ontology_term + interPro_go_ids_term))}

	return output_dic



def omicsboxPathwayExportTxt(file):
	"""
	Method to Omicsbox_pathway_export.txt format
	Args:
	     file: file path
	Returns:
			{gene: {"note": {"Pathways": value}}
 
		Not matching retuns ''
	"""
	output_dic = {}
	table = pd.read_csv(file, sep='\t', keep_default_na=False)

	for i in range(len(table)):

		pathways_term = table.iloc[i, 8]	
		print(pathways_term)

		output_dic[table.iloc[i, 0]] = {"note": {"Pathways": pathways_term}}

	return output_dic




def proteinsFaIprscnTsv(file):
	"""
	Method to Trichoderma_atroviride_IMI206040.proteins.fa.iprscn.tsv format
	Args:
	     file: file path
	Returns: 
		 Not matching retuns ['']
	"""




def xxx(file):
	"""
	Method to 
	Args:
	     file: file path
	Returns: 
		 Not matching retuns ['-']
	"""
	

