#!/usr/bin/env python3.5

entrada=open("TatrovirideIMI206040.ann.gff3","r")

lineas=entrada.readlines()

salida=open("atrovirideIMI206040.ann.fasta","w")

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
			salida.write(">"+dic_aux.get("ID","None"))
			salida.write("\n")
			salida.write(dic_aux.get("translation","None"))


			
salida.close()
entrada.close()				
		
