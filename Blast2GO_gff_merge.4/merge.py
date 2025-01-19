from filesParse import *
from goOboParse import goObo
from goCounts import goCounts

# Parse all inputs files using fileParse methods
# Gather together in one big dictionary
def merge_inputs_dictionaries(list_input_dictionaries):
	all_list_input_dictionaries={}
	for dictionary in list_input_dictionaries:
		for gene in dictionary:
			if gene in all_list_input_dictionaries.keys():
				all_list_input_dictionaries[gene]=dictionary[gene] + all_list_input_dictionaries[gene]
			else:		
				all_list_input_dictionaries[gene]=dictionary[gene]
	return all_list_input_dictionaries


def A(file_output):

	all_input_dictionaries = merge_inputs_dictionaries([
		b2GoTable("blast2go_omicsbox_table.txt", 2, 9),
		b2GoGffExport("blast2go_gff_export.gff"),
		b2GoTable("eggnog_trichoderma_atroviride_imi206040_v3_proteins_fa_eggnog_go_slim.txt", 0, 2),
		b2GoTable("eggnog_Trichoderma_atroviride_IMI206040_V3.proteins.fa_eggnog.txt", 1, 9) ,
		b2GOIprcscn("iprscan_Tatro_V3_annotations.txt"),
		b2GoAnnot("blast2go.annot")
		])

	salida = open(file_output,"w")
	salida.write("geneid"+"\t"+"GOid"+"\t"+"GOname"+"\t"+"GOnamespace"+"\n")

	dict_GO_values={}
	dic_GO_obo_values={}

	for gene in all_input_dictionaries:
		
		# Remove GO duplicated go terms
		GO_values=list(set(all_input_dictionaries[gene]))
		
		dict_GO_values[gene]=GO_values
		
		# Search for the GO termn in go.obo and retrieve GO Name and GO namespace
		# Update deprecate GO terms in the process
		GO_obo_values=goObo(GO_values)

		dic_GO_obo_values[gene]=GO_obo_values[0]

		salida.write(gene+"\t")
		salida.write(";".join(GO_obo_values[0])+"\t")
		salida.write(";".join(GO_obo_values[1])+"\t")
		salida.write(";".join(GO_obo_values[2])+"\n")

	# Los escribo en un archivo de salida
	salida.close()

	print("Output Result GO terms before obo process")
	goCounts(dict_GO_values)
	print("Output Result GO terms after obo process")
	goCounts(dic_GO_obo_values)

	GO=[]
	for key in dict_GO_values:
		GO= GO + dict_GO_values[key]
	GO = set(GO)
	GO_obo=[]
	for key in dic_GO_obo_values:
		GO_obo= GO_obo + dic_GO_obo_values[key]
	GO_obo = set(GO_obo)

	print("Difference of deprecates", len(GO - GO_obo))
	print("List of deprecates", GO - GO_obo)

A("result_go_table.txt")