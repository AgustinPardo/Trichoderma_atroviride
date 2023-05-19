#!/usr/bin/env python3.5

entrada=open("Trire_Chr_GeneCatalog_20190325.gff3","r")
lineas=entrada.readlines()
entrada.close()

salida=open("output.gff3","w")

for linea in lineas:
	features=linea.split("\t")	
	if linea[0]=="#":
		salida.write(linea)
	else:
		col_nine=features[8]
		col_nine_splitted=col_nine.split(";")	

		if "ID" in col_nine and "Name" in col_nine:
			col_nine_splitted= col_nine_splitted[1:]

		if "ID" in col_nine and "Name" not in col_nine:
			id_selection=col_nine_splitted[0].split("=")
			id_selection_value=id_selection[1]
			col_nine_splitted[0]="Name="+id_selection_value

		col_nine=";".join(col_nine_splitted)

		features[8]=col_nine
		linea="\t".join(features)
		salida.write(linea)

salida.close()