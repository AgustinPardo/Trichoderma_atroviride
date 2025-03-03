from BCBio import GFF
from filesParse import *

def print_variable_name_vars(variable):
    for name, value in vars().items():
        if value is variable:
            print(name)

"""
Example of qualifiers to be added to the GFF:
		'Tatro_000017-T1': {
			'Parent': 'Tatro_000017', 
			'product': 'hypothetical protein', 
			'Dbxref': {
				'PFAM': ['PF08659', 'PF00106', 'PF13561'], 
				'InterPro': ['IPR036291', 'IPR002347'],
                'GENE3D': ['G3DSA:3.40.50.720'], 
				'PANTHER': ['PTHR42760'], 
				'FUNFAM': ['G3DSA:3.40.50.720:FF:000084']}, 
			'note': {
				'SMCOG': ['SMCOG1001'],
				'CAZyme': ['GH93'],
				'MEROPS': ['MER0044357'],
				'BUSCO': [],
                'EC': [],
				'KEGG': ['K21248'],
                'Pathways': ['R-BTA-9845576']},
            'Ontology_term': ['GO:0005515']
            }	
"""

def add_qualifiers_to_gff(gff_input_file, 
                          additional_qualifiers, 
                          additional_qualifiers_source, 
                          gff_output_file):
    
    in_handle = open(gff_input_file)

    rec_output=[]

    for rec in GFF.parse(in_handle):
        for gene_feature in rec.features:

            for sub_features in gene_feature.sub_features:    

                if sub_features.type == "mRNA":            
                    input_gff_gene_id = sub_features.qualifiers["ID"][0]    
                    try:
                        gene_new_qualifier_dic = additional_qualifiers[input_gff_gene_id]                       
                    except:                        
                        text_gene_non_exist = input_gff_gene_id + " not exist in " + additional_qualifiers_source
                        print(text_gene_non_exist)
                        continue
                    # add the new qualifiers.

                    for new_qualifier, new_value_qualifier in gene_new_qualifier_dic.items():   
                        if new_qualifier in sub_features.qualifiers:
                            if new_qualifier in ["Dbxref", "note"]:
                                for sub_qual_key, sub_qual_value in new_value_qualifier.items():
                                    if len(sub_qual_value)>1:
                                        for element in qualifier_dic_list_format_to_gff(sub_qual_key, sub_qual_value):
                                            if element not in sub_features.qualifiers[new_qualifier]:
                                                sub_features.qualifiers[new_qualifier].append(element) 
         
                            else:   
                                if new_qualifier == "Ontology_term":
                                    for element in new_value_qualifier:
                                        if element not in sub_features.qualifiers[new_qualifier]:
                                            sub_features.qualifiers[new_qualifier].append(element) 
                                else:                    
                                    sub_features.qualifiers[new_qualifier].append(new_value_qualifier)

                        else:
                            if new_qualifier in ["Dbxref", "note"]:
                                for sub_qual_key, sub_qual_value in new_value_qualifier.items():
                                    if len(sub_qual_value)>1:
                                        sub_features.qualifiers[new_qualifier] = qualifier_dic_list_format_to_gff(sub_qual_key, sub_qual_value) 
                            else:    
                                sub_features.qualifiers[new_qualifier] = new_value_qualifier

        rec_output.append(rec)

    with open(gff_output_file, "w") as out_handle:
        GFF.write(rec_output, out_handle)

def qualifier_dic_list_format_to_gff(sub_q_key, sub_q_value):
    list_output=[]
    for element in sub_q_value:
        list_output.append(sub_q_key + ":" + element)
    return list_output

inputs_list = [
    (eggnogTxt("inputs/Trichoderma_atroviride_IMI206040_V3.proteins.fa_eggnog.txt"), "eggnogTxt_input"),
    (omicsboxTableTxt("inputs/omicsbox_table.txt"), "omicsboxTableTxt_input"),
    (omicsboxTableDiamondBlastTxt("inputs/omicsbox_table_DiamondBlast.txt"), "omicsboxTableDiamondBlastTxt_input"),
    (omicsboxPathwayExportTxt("inputs/Omicsbox_pathway_export.txt"), "omicsboxPathwayExportTxt_input"),
    (proteinsFaIprscnTsv("inputs/Trichoderma_atroviride_IMI206040.proteins.fa.iprscn.tsv"), "proteinsFaIprscnTsv_input"),
    (allAnnotationsTabular("inputs/Tatro_V3_annot_allannotations.tabular"),"allAnnotationsTabular_input"),
    (annotTbl("inputs/Tatro_V3_annot.tbl"), "annotTbl_input"),
    (antismashGbk("inputs/Tatroviride_IMI206040_antismashiV3.gbk"), "antismashGbk_input")
    ]

add_qualifiers_to_gff("inputs/Trichoderma_atroviride_IMI206040V3.gff3",
                       b2GoAnnot("inputs/blast2go.annot"),
                       "b2GoAnnot_input",
                       "outputs/enriched_output.gff")

for input in inputs_list:    
    add_qualifiers_to_gff("outputs/enriched_output.gff",
                        input[0],
                        input[1],
                        "outputs/enriched_output.gff")

#Remove remark from gff header
gff_file_lines_cleaned=[]
with open("outputs/enriched_output.gff", "r") as f:
    gff_file_lines = f.readlines()
    for line in gff_file_lines:
        if "remark" not in line:
            gff_file_lines_cleaned.append(line)
f.close()
with open("outputs/enriched_output.gff", "w") as f:
    for line in gff_file_lines_cleaned:
        if "remark" not in line:
            f.write(line)
f.close()