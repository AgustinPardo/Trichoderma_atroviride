#!python 2.7

from __future__ import print_function
from BCBio import GFF
import pandas as pd
from goatools.obo_parser import GODag

def B2Go_dic(file,reduct=True):
	table=pd.read_csv(file,sep='\t',keep_default_na=False)
	if reduct:
		table["SeqName"] = table["SeqName"].apply(lambda x : x[:15])
	dic_table=table.set_index('SeqName').T.to_dict('dict')
	return dic_table

obodag=GODag("/home/apardo/workspace/Karina/GO_merge/go.obo", optional_attrs={'consider', 'replaced_by'}, load_obsolete=True, prt=None) 

def goOboNames(go_id):
	go_name = obodag[go_id].name
	go_namespace = obodag[go_id].namespace
	return ['...'.join([go_id,go_name,go_namespace])]

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

dic_print_obsoletes={}
def goObo(lista_go_id):
	go_list_depr=[]
	go_list=[]
	for go_id in lista_go_id:
		if obodag[go_id].is_obsolete:
			go_list_depr=goOboDeprecated(go_id)+go_list_depr

			dic_print_obsoletes[go_id]=goOboDeprecated(go_id)
		else:
			go_list=[go_id]+go_list

	go_list_ND = list(dict.fromkeys(go_list+go_list_depr))

	go_list_ND_names=[]

	for go in go_list_ND:
		if type(go) is tuple:
			go_list_ND_names=go_list_ND_names+goOboNames(go[0])
		else:
			go_list_ND_names=go_list_ND_names+goOboNames(go)
	return go_list_ND_names



b2go_dic=B2Go_dic("blast2go_go_table.txt")
b2goV2_dic=B2Go_dic("blast2go_go_tableV2.txt",False)

nuevos_scaffolds = []

# Calcular
#Cuantos GO tenia el gff antes?
gff_go_before=[]
#Cuantos genes tenian al menos un GO
gff_go_gene_before=0
#Cuantos GO tiene el gff nuevo?
gff_go_after=[]
#Cuantos genes tenian al menos un GO en el nuevo gff
gff_go_gene_after=0

for scaffold in GFF.parse("scaffolds.gff"):
	for gene in scaffold.features:
		for mRNA in gene.sub_features:
			qualifier=mRNA.qualifiers

			id_gff=qualifier.get("ID",None)[0]
			hits_ref_gff=qualifier.get("hits_ref",[None])[0]
			GO_gff=qualifier.get("GO",[])


			try:
				GO_b2go = b2go_dic[id_gff]["GO IDs"].split("; ")
				GO_b2go = [x[2:] for x in GO_b2go]
			except:
				GO_b2go =[]

			try:
				#Cambio el formato a GO:###### separados por coma
				GO_b2goV2 = (b2goV2_dic[int(hits_ref_gff)]["GO IDs"].split("; "))
				GO_b2goV2 = [x[2:] for x in GO_b2goV2]
				#Cambio el formato a GO:###### separados por coma

			except:
				GO_b2goV2 =[]

			try:
				if b2goV2_dic[int(hits_ref_gff)]["InterPro GO IDs"] == "no IPS match":
					raise Exception()

				GO_InterPro_b2goV2 = (b2goV2_dic[int(hits_ref_gff)]["InterPro GO IDs"].split("; "))
				GO_InterPro_b2goV2 = [x[2:] for x in GO_InterPro_b2goV2]
				#Cambio el formato a GO:###### separados por coma

			except:
				GO_InterPro_b2goV2 = []


			GO_all = []
			GO_all = GO_gff + GO_b2go + GO_b2goV2 + GO_InterPro_b2goV2

			try:
				GO_all.remove('')
			except:
				pass
			try:
				GO_all.remove(" GO terms")
			except:
				pass

			#######
			# Cuentas
			if len(GO_gff)!=0:
				gff_go_gene_before=gff_go_gene_before+1
			if len(GO_all)!=0:
				gff_go_gene_after=gff_go_gene_after+1
	
			gff_go_before= gff_go_before + (list(dict.fromkeys(GO_gff)))
			gff_go_after= gff_go_after + (list(dict.fromkeys(GO_all)))
			#######
			
			# Saco los terminos GO repetidos
			GO_all_noDuplicates = list(dict.fromkeys(GO_all))

			if len(GO_all_noDuplicates) != 0:
				go_list=goObo(GO_all_noDuplicates)
				qualifier["GO"]=go_list


	nuevos_scaffolds.append(scaffold)

# Cuentas sobre gff antes y despues
print(len((gff_go_before)))
print(len(list(dict.fromkeys(gff_go_before))))
print(len((gff_go_after)))
print(len(list(dict.fromkeys(gff_go_after))))
print(gff_go_gene_before)
print(gff_go_gene_after)

# Rastreo de los obsoletos
# for key in dic_print_obsoletes:
# 	print (key, end=" ")
# 	print ("Consider", end=" ")
# 	if type(dic_print_obsoletes[key][0]) is tuple:
# 		print (dic_print_obsoletes[key][0])
# 	else:
# 		print (" ".join(dic_print_obsoletes[key]))


# Grabar los cambios
# with open("annotation_v2.gff","w") as h:
#     GFF.write(nuevos_scaffolds,h)

