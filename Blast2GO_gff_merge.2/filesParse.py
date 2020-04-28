#!python 2.7

import pandas as pd
from BCBio import GFF
from goatools.obo_parser import GODag

def b2GoDic(file,reduct=True):
	dic_table={}
	table=pd.read_csv(file,sep='\t',keep_default_na=False)
	if reduct:
		table["SeqName"] = table["SeqName"].apply(lambda x : x[:15])
	for i in range(len(table)):
		go_terms = table.iloc[i,9].split('; ')
		go_terms = [i[2:] for i in go_terms]
		dic_table[table.iloc[i,2]] = go_terms
	return dic_table

def iprcscnDic(file):
	dic_table={}
	table=pd.read_csv(file,sep='\t',keep_default_na=False)
	#keep_default_na=False, pone '' en vez de nan
	for i in range(len(table)):
		go_terms=table.iloc[i,13].split("|")
		if table.iloc[i,0] in dic_table.keys():
			dic_table[table.iloc[i,0]]=go_terms+dic_table[table.iloc[i,0]]
		else:
			dic_table[table.iloc[i,0]]=go_terms
	return dic_table

def goGffDic(file):
	out_dic={}
	for scaffold in GFF.parse(file):
		for gene in scaffold.features:
			for mRNA in gene.sub_features:
				qualifier=mRNA.qualifiers

				id_gff=qualifier.get("ID",None)[0]
				GO_gff=qualifier.get("Ontology_term",[''])
				out_dic[id_gff]=GO_gff
	return out_dic


obodag=GODag("/home/agustin/workspace/Karina/GO_merge/go.obo", optional_attrs={'consider', 'replaced_by'}, load_obsolete=True, prt=None) 

def goOboNames(go_id):
	go_name = obodag[go_id].name
	go_namespace = obodag[go_id].namespace
	return [[go_id],[go_name],[go_namespace]]

def goOboDeprecated(go_id):
	consider=obodag[go_id].consider
	if len(consider) != 0:
		consider_list=[]
		for term in consider:
			if term[:2] != "CL":
				consider_list=consider_list + [term]
			else:
				consider_list=consider_list + [(go_id,term)]
	else:
		# Pongo o saco los que no tiene consider
		# consider_list=[]
		consider_list=[go_id]

	return consider_list

def goObo(lista_go_id):
	go_list_depr=[]
	go_list=[]
	for go_id in lista_go_id:
		if obodag[go_id].is_obsolete:
			go_list_depr=goOboDeprecated(go_id)+go_list_depr
		else:
			go_list=[go_id]+go_list

	go_list_ND = list(dict.fromkeys(go_list+go_list_depr))

	go_list_id=[]
	go_list_name=[]
	go_list_namespace=[]
	for go in go_list_ND:
		if type(go) is tuple:
			go_list_id = goOboNames(go[0])[0]+go_list_id
			go_list_name = goOboNames(go[0])[1]+go_list_name
			go_list_namespace = goOboNames(go[0])[2]+go_list_namespace
		else:
			go_list_id = goOboNames(go)[0]+go_list_id
			go_list_name = goOboNames(go)[1]+go_list_name
			go_list_namespace = goOboNames(go)[2]+go_list_namespace

	return 	go_list_id, go_list_name, go_list_namespace
