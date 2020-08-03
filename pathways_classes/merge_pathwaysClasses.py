
#Leo archivos de macro clases y creo una lista
def macroClass(file):
    entrada=open(file)
    lineas=entrada.readlines()
    lista_pathways=[]
    for linea in lineas:
        if linea!="\n":
            lista_pathways.append(linea.split("\t")[0])
    entrada.close()
    return lista_pathways

biosynthesis_list=macroClass("/home/agustin/workspace/Karina/pathways_classes/biosynthesis.txt")


def microClass(file):
    entrada=open(file)
    lineas=entrada.readlines()
    dic_pathways={}
    for linea in lineas:
        if linea!="\n":
            dic_pathways[(linea.split("\t")[0])]=linea.split("\t")[2].rstrip()
    entrada.close()
    return dic_pathways

# Leo el archivo de genes y coordenadas
def genesCoordinates(file):
    entrada=open(file)
    lineas=entrada.readlines()
    dic_genes={}
    for linea in lineas:
        if linea!="\n":
            dic_genes[(linea.split("\t")[0])]=linea.split("\t")[1:]
    entrada.close()
    return dic_genes

microClasses=microClass("/home/agustin/workspace/Karina/pathways_classes/class-pathways-2.txt")

activation_inactivation_list=macroClass("/home/agustin/workspace/Karina/pathways_classes/activation-inactivation.txt")
biosynthesis_list=macroClass("/home/agustin/workspace/Karina/pathways_classes/biosynthesis.txt")
degradation_list=macroClass("/home/agustin/workspace/Karina/pathways_classes/Degradation.txt")
detoxification_list=macroClass("/home/agustin/workspace/Karina/pathways_classes/Detoxification.txt")
energy_metabolism_list=macroClass("/home/agustin/workspace/Karina/pathways_classes/Energy-metabolism.txt")
glycan_pathways_list=macroClass("/home/agustin/workspace/Karina/pathways_classes/Glycan-pathways.txt")
macromolecule_modification_list=macroClass("/home/agustin/workspace/Karina/pathways_classes/Macromolecule-Modification.txt")
metabolic_clusters_list=macroClass("/home/agustin/workspace/Karina/pathways_classes/Metabolic-clusters.txt")

dic_genes_coordinate=genesCoordinates("/home/agustin/workspace/Karina/pathways_classes/allgenes-pathways.txt")

entrada=open("/home/agustin/workspace/Karina/pathways_classes/pathways-3.txt","r")
salida=open("micro_macro_coordinates_pathways.txt","w")
lineas=entrada.readlines()
salida.write("frame"+"\t"+"name"+"\t"+"microClass"+"\t"+"macroClass"+"\t"+"genes"+"\t"+"coordinates"+"\n")
for linea in lineas:
    if linea!="\n" and linea[0:5]!= "Frame":

        linea_split=linea.split("\t")

        pathway_id=linea_split[0]
        pathway_name=linea_split[1]
        pathway_genes=linea_split[2]
        pathway_microClass=microClasses[pathway_id]

        salida.write(pathway_id+"\t"+pathway_name+"\t"+pathway_microClass+"\t")

        macroClass_list=[]

        if pathway_id in activation_inactivation_list:
            macroClass_list.append("activation-inactivation")
        if pathway_id in biosynthesis_list:
            macroClass_list.append("biosynthesis")
        if pathway_id in degradation_list:
            macroClass_list.append("degradation")
        if pathway_id in detoxification_list:
            macroClass_list.append("detoxification")
        if pathway_id in energy_metabolism_list:
            macroClass_list.append("energy-metabolism")
        if pathway_id in glycan_pathways_list:
            macroClass_list.append("glycan-pathways")
        if pathway_id in macromolecule_modification_list:
            macroClass_list.append("macromolecule-modification")
        if pathway_id in metabolic_clusters_list:
            macroClass_list.append("metabolic-clusters")

        salida.write(";".join(macroClass_list)+"\t")

        pathway_genes_list=pathway_genes.rstrip("\"").rstrip("\n").rstrip("\"").lstrip("\"").split(", ")
        salida.write(";".join(pathway_genes_list)+"\t")

        coordinates_list=[]
        for gene in pathway_genes_list:
            if gene == "":
                pass
            if gene[0:5] == "G23KX":
                coordinates_list.append("None"+" "+"None")
            else:
                coordinates=dic_genes_coordinate.get(gene,"None")
                coordinates_list.append(coordinates[0]+" "+coordinates[1])
        salida.write(";".join(coordinates_list))
        salida.write("\n")

salida.close()
entrada.close()
