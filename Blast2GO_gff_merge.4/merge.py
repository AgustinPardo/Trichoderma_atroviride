from filesParse import *
from goOboParse import goObo


def merge_dictionaries(list_input_dictionaries):
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
def merge_input_dictionaries(*args):
	# Parse all inputs files using fileParse methods
	# Gather together in one big dictionary	
	return merge_dictionaries([
		b2GoTable("inputs/blast2go_omicsbox_table.txt", 2, 9),
		b2GoGffExport("inputs/blast2go_gff_export.gff"),
		b2GoTable("inputs/eggnog_trichoderma_atroviride_imi206040_v3_proteins_fa_eggnog_go_slim.txt", 0, 2),
		b2GoTable("inputs/eggnog_Trichoderma_atroviride_IMI206040_V3.proteins.fa_eggnog.txt", 1, 9),
		b2GOIprcscn("inputs/iprscan_Tatro_V3_annotations.txt"),
		b2GoAnnot("inputs/blast2go.annot")
		])


def get_set_GO_from_dict_gene_GO(dict_gene_GO):
	GO=[]
	for key in dict_gene_GO:
		GO = GO + dict_gene_GO[key]
	return set(GO)


def compare_GO_from_dict_gene_GO(dict1, dict2):
	GO=get_set_GO_from_dict_gene_GO(dict1)
	GO_obo = get_set_GO_from_dict_gene_GO(dict2)
	print("Difference of deprecates", len(GO - GO_obo))
	print("List of deprecates", GO - GO_obo)


@outputReview
def B(*args):

	all_input_dictionaries = merge_input_dictionaries("Output Result GO terms before obo process")

	# Store results in a file
	salida = open("outputs/results_go_table.txt","w")
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

	compare_GO_from_dict_gene_GO(all_input_dictionaries, GO_values_obo_counts)

	return GO_values_obo_counts

GO_values_obo_counts = B("Output Result GO terms after obo process")