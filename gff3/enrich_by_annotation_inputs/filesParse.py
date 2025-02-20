import pandas as pd
import re

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
	Method to Trichoderma_atroviride_IMI206040_V3.proteins.fa_eggnog.txt format
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
		{'Tatro_005317': {'Pathways': [''], 'Ontology_term': [''], 'Dbxref': {'Gene3D': ['G3DSA:3.40.50.720'], 'PANTHER': ['PTHR42760'], 'Pfam': ['PF13561'], 'FunFam': ['G3DSA:3.40.50.720:FF:000084']}}}
		 Not matching retuns ['']
	"""
	output_dic = {}

	table = pd.read_csv(file, sep='\t', keep_default_na=False)

	for i in range(len(table)):
		
		gene = table.iloc[i, 0]

		dbxref_term_type = table.iloc[i, 3]
		dbxref_term = table.iloc[i, 4]

		if dbxref_term_type in ["PANTHER", "FunFam", "Pfam","ProSiteProfiles","Gene3D"]:
			if gene in output_dic.keys():
				if "Dbxref" in output_dic[gene].keys():
					if dbxref_term_type in output_dic[gene]["Dbxref"]:
						output_dic[gene]["Dbxref"].update({dbxref_term_type: list(set(output_dic[gene]["Dbxref"][dbxref_term_type] + [dbxref_term]))})
					else:
						output_dic[gene]["Dbxref"].update({dbxref_term_type: [dbxref_term]})
				else:
					output_dic[gene].update({"Dbxref": {dbxref_term_type: [dbxref_term]}})
			else:
				output_dic.update({gene: {"Dbxref": {dbxref_term_type: [dbxref_term]}}})


		pathway_term = table.iloc[i, 14]
		pathway_term = pathway_term.replace("-", "")
		if gene in output_dic.keys():
			if "Pathways" in output_dic[gene]:
				output_dic[gene]["Pathways"] = list(set(output_dic[gene]["Pathways"] + [pathway_term]))
			else:
				output_dic[gene]["Pathways"] = [pathway_term]
		else:
			output_dic[gene] = {"Pathways": [pathway_term]}


		ontology_term = table.iloc[i, 13]
		ontology_term = ontology_term.replace("-", "")
		if gene in output_dic.keys():
			if "Ontology_term" in output_dic[gene]:
				output_dic[gene]["Ontology_term"] = list(set(output_dic[gene]["Ontology_term"] + [ontology_term]))
			else:
				output_dic[gene]["Ontology_term"] = [ontology_term]
		else:
			output_dic[gene] = {"Ontology_term": [ontology_term]}	

	return output_dic


def allAnnotationsTabular(file):
	"""
	Method to Tatro_V3_annot_allannotations.tabular format
	Args:
	     file: file path
	Returns: 
		 Not matching retuns ['-']
	"""

	output_dic = {}

	table = pd.read_csv(file, sep='\t', keep_default_na=False)

	for i in range(len(table)):
		
		gene = table.iloc[i, 1]
		parent_term = table.iloc[i, 0]
		product_term = table.iloc[i, 7]
		dbxref_pfam_term = table.iloc[i,12]

		dbxref_interPro_term = table.iloc[i,13]
		pattern = r"\bIPR\w+\b"
		dbxref_interPro_term = re.findall(pattern, dbxref_interPro_term)

		note_ECnumber_term = table.iloc[i,10]
		note_AntiSMASH_term = table.iloc[i,21]

		note_notes_term = table.iloc[i,22]
		pattern = r"\bSMCOG\w+\b"
		note_notes_term = re.findall(pattern, note_notes_term)

		ontology_term = table.iloc[i,16]
		print(ontology_term)
		# print(gene, parent_term, product_term, dbxref_pfam_term, 
		# dbxref_interPro_term, note_ECnumber_term, note_AntiSMASH_term, note_notes_term,
		# ontology_term, sep="*")
	
	return output_dic

allAnnotationsTabular("inputs/Tatro_V3_annot_allannotations.tabular")


def xxxxx(file):
	"""
	Method to 
	Args:
	     file: file path
	Returns: 
		 Not matching retuns ['-']
	"""
	pass
	
