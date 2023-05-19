#!/usr/bin/env python3.5

entrada=open("/home/apardo/workspace/Karina/Blast2GO_gff_merge/annotation_v2.gff","r")

lineas=entrada.readlines()

salida=open("ID_GOes.txt","w")

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
				
			salida.write(dic_aux.get("ID","None")+"\t")
			salida.write(dic_aux.get("hits_ref","None")+"\t")

			GOes=dic_aux.get("GO","None")

			GO_term =[]
			GO_name =[]
			GO_namespace =[]

			if GOes == "None":
				GO_term="None"
				GO_name="None"
				GO_namespace="None"
				salida.write(GO_term+"\t"+GO_name+"\t"+GO_namespace)
			else:
				GOes=GOes.split(",")
				for GO in GOes:
					GO_elemt=GO.split("...")
					GO_term.append(GO_elemt[0]) 
					GO_name.append(GO_elemt[1])
					GO_namespace.append(GO_elemt[2])
				salida.write(",".join(GO_term)+"\t"+",".join(GO_name)+"\t"+",".join(GO_namespace))

			salida.write("\n")
			
salida.close()
entrada.close()				
		
