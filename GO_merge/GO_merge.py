#!/usr/bin/env Python 3.7.0


from goatools.obo_parser import GODag

locus_tags_file=open("locus_tags_data.txt","r")
locus_tags_lines= locus_tags_file.readlines()

obodag = GODag("go.obo")

salida=open("locus_tags_GO_merge.txt","w")	


for line in locus_tags_lines:
	line_split=line.split("\t")
	GO_term=line_split[1]
	GO_terms=GO_term.split(",")
	
	salida.write(line_split[0]+"\t")
	salida.write(line_split[2]+"\t")
	salida.write(line_split[1]+"\t")

	name_list=[]
	namespace_list=[]
	for term in GO_terms:
		try:
			name_list.append(obodag.get(term).name)
		except:
			name_list.append("None")
	salida.write(",".join(name_list))
	
	salida.write("\t")
	
	for term in GO_terms:
		try:
			namespace_list.append(obodag.get(term).namespace)
		except:
			namespace_list.append("None")
	salida.write(",".join(namespace_list))
	salida.write("\n")
	
locus_tags_file.close()
salida.close()	
	

