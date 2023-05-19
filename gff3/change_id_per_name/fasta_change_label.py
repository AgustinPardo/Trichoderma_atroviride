#!/usr/bin/env python3.5

entrada=open("Trire_Chr_GeneCatalog_CDS_20190325.fasta","r")
lineas=entrada.readlines()
entrada.close()
salida=open("output.fasta","w")

# FROM "jgi|Trire_Chr|109159|TrA0001W" to "jgi.p|Trire_Chr|108761"
for linea in lineas:
	if linea[0]==">":
		linea=linea[1:].rstrip("\n")
		linea_split=linea.split("|")
		linea="|".join([">jgi.p","Trire_Chr",linea_split[2]])+"\n"
		salida.write(linea)
	else:
		salida.write(linea)

salida.close()		