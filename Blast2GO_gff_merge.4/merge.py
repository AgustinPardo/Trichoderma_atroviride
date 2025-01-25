from filesParse import *
from goOboParse import goObo


def merge_inputs_dictionaries(list_input_dictionaries):
	all_list_input_dictionaries={}
	for dictionary in list_input_dictionaries:
		for gene in dictionary:
			if gene in all_list_input_dictionaries.keys():
				all_list_input_dictionaries[gene]=dictionary[gene] + all_list_input_dictionaries[gene]
				# Remove duplicated GO terms
				all_list_input_dictionaries[gene] = list(set(all_list_input_dictionaries[gene]))
			else:		
				all_list_input_dictionaries[gene]=dictionary[gene]
	return all_list_input_dictionaries

@outputReview
def A(*args):
	# Parse all inputs files using fileParse methods
	# Gather together in one big dictionary	
	return merge_inputs_dictionaries([
		b2GoTable("blast2go_omicsbox_table.txt", 2, 9),
		b2GoGffExport("blast2go_gff_export.gff"),
		b2GoTable("eggnog_trichoderma_atroviride_imi206040_v3_proteins_fa_eggnog_go_slim.txt", 0, 2),
		b2GoTable("eggnog_Trichoderma_atroviride_IMI206040_V3.proteins.fa_eggnog.txt", 1, 9),
		b2GOIprcscn("iprscan_Tatro_V3_annotations.txt"),
		b2GoAnnot("blast2go.annot")
		])

@outputReview
def B(file_output):

	all_input_dictionaries = A("Output Result GO terms before obo process")

	# Store results in a file
	salida = open(file_output,"w")
	salida.write("geneid"+"\t"+"GOid"+"\t"+"GOname"+"\t"+"GOnamespace"+"\n")

	GO_values_obo_counts={}

	for gene in all_input_dictionaries:
		
		GO_values=all_input_dictionaries[gene]
				
		# Search for the GO termn in go.obo and retrieve GO Name and GO namespace
		# Update deprecate GO terms in the process
		GO_values_obo=goObo(GO_values)

		salida.write(gene+"\t")
		salida.write(";".join(GO_values_obo[0])+"\t")
		salida.write(";".join(GO_values_obo[1])+"\t")
		salida.write(";".join(GO_values_obo[2])+"\n")
		
		GO_values_obo_counts[gene]=GO_values_obo[0]

	salida.close()

	# GO=[]
	# for key in dict_GO_values:
	# 	GO= GO + dict_GO_values[key]
	# GO = set(GO)
	# GO_obo=[]
	# for key in dic_GO_obo_values:
	# 	GO_obo= GO_obo + dic_GO_obo_values[key]
	# GO_obo = set(GO_obo)

	# print("Difference of deprecates", len(GO - GO_obo))
	# print("List of deprecates", GO - GO_obo)

	return GO_values_obo_counts

B("result_go_table.txt")