from Bio import SeqIO
import pandas as pd
import re

def b2GoAnnot(file):
	"""
	Method to blast2go.annot format
	Args:
	     file: file path
	Returns: 
	     'Tatro_012421-T1': {
			'Ontology_term': ['GO:0003678', 'GO:0000723', 'GO:0006281']}

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
			output_dic[gene] = {"Ontology_term": [ontology_term]}

	return output_dic

# print(b2GoAnnot("inputs/blast2go.annot"))

def eggnogTxt(file):
	"""
	Method to Trichoderma_atroviride_IMI206040_V3.proteins.fa_eggnog.txt format
	Args:
	     file: file path
	Returns: 
	    'Tatro_012437-T1': {
			'product': 'endoplasmic reticulum organization',
			'note': {
				'EC': [],
				'KEGG': ['K21248']},
			'Ontology_term': ['GO:0006971', 'GO:0031152']}

		Not matching retuns []
	"""
	output_dic = {}
	table = pd.read_csv(file, sep='\t', keep_default_na=False)

	for i in range(len(table)):
		product_term = table.iloc[i, 3]	

		note_EC_term = [table.iloc[i, 7]]
		note_EC_term.remove("") if "" in note_EC_term else None

		note_KEGG_term = [table.iloc[i, 11]]
		note_KEGG_term.remove("") if "" in note_KEGG_term else None

		ontology_term = table.iloc[i, 9].split('; ')
		ontology_term = [i[2:] for i in ontology_term]

		if "" in ontology_term:
			ontology_term.remove("")

		output_dic[table.iloc[i, 1]] = {
		"product": product_term,
		"note": {"EC": note_EC_term, "KEGG": note_KEGG_term},
		"Ontology_term": ontology_term}

	return output_dic

# print(eggnogTxt("inputs/Trichoderma_atroviride_IMI206040_V3.proteins.fa_eggnog.txt"))

def omicsboxTableTxt(file):
	"""
	Method to omicsbox_table.txt format
	Args:
	     file: file path
	Returns: 
	    'Tatro_012399-T1': {
			'Dbxref': {
				'InterPro': ['IPR001138', 'IPR001138']},
				'PFAM': ['PF00106', 'PF00107']} 
			'Ontology_term': ['GO:0016787', 'GO:0005737']}

		Not matching retuns []
	"""
	output_dic = {}
	table = pd.read_csv(file, sep='\t', keep_default_na=False)

	for i in range(len(table)):
		dbxref_term = table.iloc[i, 13]	
		dbxref_term = dbxref_term.split(";")
		for term in dbxref_term:
			if term not in ["no IPS match"]:
				dbxref_term = table.iloc[i,13]
				pattern_IPR = r"\bIPR\w+\b"
				pattern_PF = r"\bPF\w+\b"
				dbxref_interPro_term = re.findall(pattern_IPR, dbxref_term)
				dbxref_pfam_term = re.findall(pattern_PF, dbxref_term)
				dbxref_pfam_term = [elem for elem in dbxref_pfam_term if elem not in ["PFAM"]]
    
		ontology_term = table.iloc[i, 9].split('; ')
		ontology_term = [i[2:] for i in ontology_term]

		ontology_term.remove("") if "" in ontology_term else None		

		interPro_go_ids_term=table.iloc[i, 14].split('; ')
		interPro_go_ids_term = [i[2:] for i in interPro_go_ids_term]

		interPro_go_ids_term.remove("no GO terms") if "no GO terms" in interPro_go_ids_term else None
		interPro_go_ids_term.remove(" GO terms") if " GO terms" in interPro_go_ids_term else None
		interPro_go_ids_term.remove(" IPS match") if " IPS match" in interPro_go_ids_term else None
		
		output_dic[table.iloc[i, 2]] = {
			"Dbxref":{"InterPro": list(set(dbxref_interPro_term)), "PFAM": list(set(dbxref_pfam_term))},
			"Ontology_term": list(set(ontology_term + interPro_go_ids_term))}
	return output_dic

# print(omicsboxTableTxt("inputs/omicsbox_table.txt"))


def omicsboxTableDiamondBlastTxt(file):
	"""
	Method to omicsbox_table_DiamondBlast.txt format
	Args:
	     file: file path
	Returns:
		'Tatro_012389-T1': {
			'Dbxref': {
				'InterPro': ['IPR001138', 'IPR001138']},
			'note': {
				'EC': []}, 
			'Ontology_term': ['GO:0005634', 'GO:0006355']}
 
		Not matching retuns ''
	"""
	output_dic = {}
	table = pd.read_csv(file, sep='\t', keep_default_na=False)

	for i in range(len(table)):
		dbxref_term = table.iloc[i, 13]	

		dbxref_term = dbxref_term.split(" ")
		interPro_term = [term for term in dbxref_term if "IPR" in term]

		note_enzyme_code_term = table.iloc[i, 11].split('; ')
		note_enzyme_code_term = [i for i in note_enzyme_code_term]
		note_enzyme_code_term.remove("") if "" in note_enzyme_code_term else None		

		ontology_term = table.iloc[i, 9].split('; ')
		ontology_term = [i[2:] for i in ontology_term]
		ontology_term.remove("") if "" in ontology_term else None		

		interPro_go_ids_term=table.iloc[i, 14].split('; ')
		interPro_go_ids_term = [i[2:] for i in interPro_go_ids_term]

		interPro_go_ids_term.remove("no GO terms") if "no GO terms" in interPro_go_ids_term else None
		interPro_go_ids_term.remove(" GO terms") if " GO terms" in interPro_go_ids_term else None
		interPro_go_ids_term.remove(" IPS match") if " IPS match" in interPro_go_ids_term else None
		
		output_dic[table.iloc[i, 2]] = {
			"Dbxref": {"InterPro": interPro_term},
			"note": {"EC": note_enzyme_code_term},
			"Ontology_term": list(set(ontology_term + interPro_go_ids_term))}

	return output_dic

# print(omicsboxTableDiamondBlastTxt("inputs/omicsbox_table_DiamondBlast.txt"))

def omicsboxPathwayExportTxt(file):
	"""
	Method to Omicsbox_pathway_export.txt format
	Args:
	     file: file path
	Returns:
		'Tatro_009201-T1': 
			{'note': {
				'Pathways': ['R-BTA-9845576']}
 
		Not matching retuns ''
	"""
	output_dic = {}
	table = pd.read_csv(file, sep='\t', keep_default_na=False)

	for i in range(len(table)):
		note_pathways_term = table.iloc[i, 8]	
		output_dic[table.iloc[i, 0]] = {
			"note": {"Pathways": [note_pathways_term]}}

	return output_dic

# print(omicsboxPathwayExportTxt("inputs/Omicsbox_pathway_export.txt"))


def proteinsFaIprscnTsv(file):
	"""
	Method to Trichoderma_atroviride_IMI206040.proteins.fa.iprscn.tsv format
	Args:
	     file: file path
	Returns: 
		Tatro_005317-T1': {
			'Dbxref': {
				'InterPro': ['IPR002347', 'IPR036291', 'IPR020904'],
				'GENE3D': ['G3DSA:3.40.50.720'], 
				'PANTHER': ['PTHR42760'], 
				'PFAM': ['PF13561'], 
				'FUNFAM': ['G3DSA:3.40.50.720:FF:000084']}}

		 Not matching retuns []
	"""
	output_dic = {}

	table = pd.read_csv(file, sep='\t', keep_default_na=False)

	for i in range(len(table)):
		
		gene = table.iloc[i, 0]

		dbxref_term_type = table.iloc[i, 3]
		dbxref_term = table.iloc[i, 4]

		if dbxref_term_type in ["PANTHER", "FunFam", "Pfam", "ProSiteProfiles", "Gene3D"]:
			dbxref_term_type=dbxref_term_type.upper()
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


		dbxref_interPro_term = table.iloc[i, 11]
		dbxref_interPro_term = dbxref_interPro_term.replace("-", "")

		if gene in output_dic.keys():
			if "Dbxref" in output_dic[gene].keys():
				if "InterPro" in output_dic[gene]["Dbxref"]:
					interPro_list = list(set(output_dic[gene]["Dbxref"]["InterPro"] + [dbxref_interPro_term]))
					if "" in interPro_list: interPro_list.remove("")
					output_dic[gene]["Dbxref"].update({"InterPro": interPro_list})
				else:
					output_dic[gene]["Dbxref"].update({"InterPro": [dbxref_interPro_term]})
			else:
				output_dic[gene].update({"Dbxref": {"InterPro": [dbxref_interPro_term]}})
		else:
			output_dic.update({gene: {"Dbxref": {"InterPro": [dbxref_interPro_term]}}})

	return output_dic

# print(proteinsFaIprscnTsv("inputs/Trichoderma_atroviride_IMI206040.proteins.fa.iprscn.tsv"))


def allAnnotationsTabular(file):
	"""
	Method to Tatro_V3_annot_allannotations.tabular format
	Args:
	     file: file path
	Returns: 
	 'Tatro_012365-T1': {
	 	'Parent': 'Tatro_012365',
		'product': 'FUB1', 
		'Dbxref': {
			'InterPro': ['IPR001227', 'IPR006162'],
			'PFAM': ['PF00106', 'PF00107']},
			'note': {'SMCOG': ['SMCOG1002'],
			  		 'BUSCO':['EOG0926390Q'],
					 'CAZyme': []}
		}
		 Not matching retuns []
	"""

	output_dic = {}

	table = pd.read_csv(file, sep='\t', keep_default_na=False)

	for i in range(len(table)):
		
		gene = table.iloc[i, 1]
		parent_term = table.iloc[i, 0]
		product_term = table.iloc[i, 7]

		dbxref_pfam_term = table.iloc[i,12]
		dbxref_pfam_term = dbxref_pfam_term.split(";")
		dbxref_pfam_term.remove("") if "" in dbxref_pfam_term else None	

		dbxref_interPro_term = table.iloc[i,13]
		pattern_IPR = r"\bIPR\w+\b"
		dbxref_interPro_term = re.findall(pattern_IPR, dbxref_interPro_term)

		busco_column = table.iloc[i,11]
		note_BUSCO_term = busco_column.split(";")
		note_BUSCO_term.remove("") if "" in note_BUSCO_term else None	

		cazyme_column = table.iloc[i,20]
		note_CAZyme_term = cazyme_column.split(";")
		note_CAZyme_term.remove("") if "" in note_CAZyme_term else None	

		notes_column = table.iloc[i,22]
		pattern_SMCOG = r"\bSMCOG\w+\b"
		note_SMCOG_term = re.findall(pattern_SMCOG, notes_column)	

		output_dic[gene] = {
				"Parent": parent_term,
				"product": product_term,
				"Dbxref":{"InterPro": dbxref_interPro_term, "PFAM": dbxref_pfam_term},
				"note": {"BUSCO": note_BUSCO_term, "CAZyme": note_CAZyme_term, "SMCOG": note_SMCOG_term,}}
	
	return output_dic

# print(allAnnotationsTabular("inputs/Tatro_V3_annot_allannotations.tabular"))

def annotTbl(file):
	"""
	Method to Tatro_V3_annot.tbl format
	Args:
	     file: file path
	Returns: 
		'Tatro_000017-T1': {
			'Parent': 'Tatro_000017', 
			'product': 'hypothetical protein', 
			'Dbxref': {
				'PFAM': ['PF08659', 'PF00106', 'PF13561'], 
				'InterPro': ['IPR036291', 'IPR002347']}, 
			'note': {
				'SMCOG': ['SMCOG1001'],
				'CAZy': ['GH93'],
				'MEROPS': ['MER0044357'],
				'BUSCO': []}}				

		 Not matching retuns []
	"""
	output_dic = {}

	file = open(file)
	lines = file.readlines()

	writer_check=False

	for row in lines:	

		row_split=row.split("\t")

		try:
			if row_split[2].rstrip() == "CDS":
				writer_check=True
				dbxref_term_dic={}
				note_term_dic={}

			if row_split[3] == "locus_tag":
				parent_term = row_split[4].rstrip()
			if row_split[3] == "protein_id":
				gene = row_split[4].split("|")[2].rstrip()
			if row_split[3] == "product":
				product_term = row_split[4].rstrip()
			
			if row_split[3] == "db_xref":
				dbxref_term = row_split[4].rstrip().split(":")
				dbxref_term_key = dbxref_term[0].rstrip()
				dbxref_term_value = dbxref_term[1].rstrip()
				if dbxref_term_key in dbxref_term_dic.keys():				
					dbxref_term_dic[dbxref_term_key] =  list(set(dbxref_term_dic[dbxref_term_key] + [dbxref_term_value]))
				else:
					dbxref_term_dic[dbxref_term_key] = [dbxref_term_value]
			
			if row_split[3] == "note":
				note_term = row_split[4].rstrip().split(":")
				note_term_key = note_term[0].rstrip()
				note_term_value = note_term[1].rstrip()

				if note_term_key != "antiSMASH":
					if "SMCOG" in note_term_key:
						note_term_value=note_term_key
						note_term_key="SMCOG"			

					if "CAZy" == note_term_key:
						note_term_key="CAZyme"

					if dbxref_term_key in note_term_dic.keys():				
						note_term_dic[note_term_key] =  list(set(note_term_dic[note_term_key] + [note_term_value]))
					else:
						note_term_dic[note_term_key] = [note_term_value]

			if writer_check:
				output_dic[gene] = {
					"Parent": parent_term,
					"product": product_term,
					"Dbxref": dbxref_term_dic,
					"note": note_term_dic}
				writer_check=False

		except IndexError:
			pass

	return output_dic
	
# annotTbl("inputs/Tatro_V3_annot.tbl")

def antismashGbk(file):
	"""
	Method to Tatroviride_IMI206040_antismashiV3.gbk format
	Args:
	     file: file path
	Returns: 
		'Tatro_010634-T1': {
			'Parent': 'Tatro_010634', 
			'product': 'hypothetical protein', 
			'Dbxref': {
				'PFAM': ['PF00498.29']},
			'Ontology_term': ['GO:0005515']}

		 Not matching retuns []
	"""
	output_dic = {}
	gbk_input = SeqIO.parse(file, "genbank") 
	for seq_record in gbk_input:

		passing_feature_checker=[]
		parent_gene_term=""
		parent_CDS_term=""

		for feature in seq_record.features:
	
			if feature.type == "gene":
				passing_feature_checker.append("gene")
				parent_gene_term = feature.qualifiers["locus_tag"][0]

			if feature.type == "CDS":
				passing_feature_checker.append("CDS")		
				gene = feature.qualifiers["protein_id"][0][5:]
				parent_CDS_term = feature.qualifiers["locus_tag"][0]
				product_term = feature.qualifiers["product"][0]

				dbxref_pfam_term=[]
				ontology_term_parsed=[]

			if feature.type == "PFAM_domain":
				dbxref_term = feature.qualifiers["db_xref"]
				dbxref_pfam_term = list(set([term for term in dbxref_term if "PF" in term]+dbxref_pfam_term))
				
				if "gene_ontologies" in feature.qualifiers.keys():
					ontology_term = feature.qualifiers["gene_ontologies"]
					ontology_term_parsed = list(set([term.split(": ")[0] for term in ontology_term]+ontology_term_parsed))	

			if (parent_gene_term!=parent_CDS_term 
				and "CDS" in passing_feature_checker 
				and "gene" in passing_feature_checker):
				output_dic[gene] = {
						"Parent": parent_CDS_term,
						"product": product_term,
						"Dbxref":{"PFAM": dbxref_pfam_term},
						"Ontology_term": ontology_term_parsed}
				passing_feature_checker=[]
	print(output_dic["Tatro_002739-T1"])
	return output_dic

# antismashGbk("inputs/Tatroviride_IMI206040_antismashiV3.gbk")

