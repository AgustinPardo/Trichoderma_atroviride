#!python 2.7

from filesParse import *

# Extraigo los terminso GO de cada archivo
# Diccionarios con una lista de terminos GO, 'Tatro_000998-T1': ['GO:0006333', 'GO:0005634']
b2go_IMI_dic=b2GoDic("Tatroviride_IMI_B2GO.txt",True)
b2go_annotation_dic=b2GoDic("TatrovirideIMI206040-GO-annotation.txt",True)
iprscan_dic=iprcscnDic("iprscan_2.tsv")
gff_dic=goGffDic("Trichoderma_atroviride_IMI206040.gff3")

dic_list=[b2go_IMI_dic,b2go_annotation_dic,iprscan_dic,gff_dic]

# Para calcular las cuentitas sobre un diccionario con listas de GO.
def cuentitas(dic_enter):
	gen_conGO=0
	gen_sinGO=0
	gen_conmasdeunGO=0
	GO_totales=[]

	for key in dic_enter:
			elements= list(set(dic_enter[key]))
			for elemento in elements:
				if elemento not in GO_totales:
					GO_totales.append(elemento)
			if len(elements) == 1 and '' in elements:
				gen_sinGO=1+gen_sinGO
			else:
				if len(elements)==1:
					gen_conGO=1+gen_conGO
				else:
					gen_conmasdeunGO=1+gen_conmasdeunGO	

	print(gen_conGO)
	print(gen_conmasdeunGO)
	print(gen_sinGO)
	print(len(dic_enter))
	print(len(GO_totales))

cuentitas(b2go_IMI_dic)

# Los junto en un solo diccionario
total_dic={}
for dictionary in dic_list:
	for key in dictionary:

		if key in total_dic.keys():
			total_dic[key]=dictionary[key]+total_dic[key]
		else:		
			total_dic[key]=dictionary[key]

# Saco los que son iguales
# Busco el GO name y GOspacename en go.obo
# Actualizo los viejos

gen_conGO=0
gen_sinGO=0
gen_conmasdeunGO=0
GO_totales=[]

salida = open("merge.txt","w")
salida.write("geneid"+"\t"+"GOid"+"\t"+"GOname"+"\t"+"GOnamespace"+"\n")
for key in total_dic:
	values=list(set(total_dic[key]))

	if '' in values:
		values.remove('')
	values=goObo(values)

	if len(values[0])==0:
		gen_sinGO=gen_sinGO+1

	if len(values[0])==1:
		gen_conGO=gen_conGO+1

	if len(values[0])>1:
		gen_conmasdeunGO=gen_conmasdeunGO+1

	for elemento in values[0]:
		if elemento not in GO_totales:
				GO_totales.append(elemento)

	salida.write(key+"\t")
	salida.write(";".join(values[0])+"\t")
	salida.write(";".join(values[1])+"\t")
	salida.write(";".join(values[2])+"\n")

# Los escribo en un archivo de salida
salida.close()

# Cuentidas del merge
print(gen_conGO)
print(gen_conmasdeunGO)
print(gen_sinGO)
print(len(total_dic))
print(len(GO_totales))

# b2go_IMI_list=list(dict.fromkeys(b2go_IMI_dic))
# b2go_annotation_list=list(dict.fromkeys(b2go_annotation_dic))
# iprscan_list=list(dict.fromkeys(iprscan_dic))
# gff_list=list(dict.fromkeys(gff_dic))

# library
# import matplotlib.pyplot as plt
# from matplotlib_venn import venn3
# venn3([set(b2go_IMI_list), set(b2go_annotation_list), set(iprscan_list), set(gff_list)],set_labels = ('b2go_IMI', 'b2go_annotation', 'iprscan', 'gff'))
# plt.show()