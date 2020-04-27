#!/usr/bin/env python3.5

import sys

#Cargo el archivo gff
ids_file=open(sys.argv[1])#open("/home/agustin/workspace/karina/gff3_extraction/TatrovirideIMI206040.ann.gff3","r")
lineas=ids_file.readlines()
dic_old_new={}

for linea in lineas:
	features=linea.split("\t")	
	if linea[0]=="#":
		pass
	else:	
		if features[2]=="mRNA":
			target=features[8].split(";")
			dic_aux={}
			for x in target:	
				x_split=x.split("=")
				dic_aux[x_split[0]]=x_split[1]
				
			dic_old_new[dic_aux.get("hits_ref","None")]=dic_aux.get("ID")
			
ids_file.close()

#Cargo el archivo de datos transcriptomicos (tsv)
ids_file=open(sys.argv[1
ed_file=open(sys.argv[2])#open("FDR2-I30-c.txt","r")
lineas=ed_file.readlines()

salida=open(sys.argv[2].split(".")[0]+"_ids_change.txt","w")

for linea in lineas:
	cols=linea.split("	")
	if linea[0]=="	":
		salida.write("\t"+"\t".join(cols[1:]))
	else:		
		print(cols[0])		
		salida.write(dic_old_new.get(cols[0],cols[0])+"\t"+ "\t".join(cols[1:]))
		
salida.close()
ids_file.close()
ed_file.close()	
		


