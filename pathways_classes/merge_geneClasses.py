entrada=open("/home/agustin/workspace/Karina/pathways_classes/micro_macro_coordinates_pathways.txt","r")
salida =open("micro_macro_coordinates_genes.txt","w")

lines=entrada.readlines()

salida.write("gene"+"\t"+"position"+"\t"+"contig"+"\t"+"frame"+"\t"+"frameName"+"\t"+"microClass"+"\t"+"macroClass"+"\n")
for line in lines[1:]:
    split_line=line.split("\t")
    genes=split_line[4].split(";")
    pathway=split_line[0]
    pathway_name=split_line[1]
    microClass=split_line[2]
    macroClass=split_line[3]
    coordinates=split_line[5].split(";")
    i=0
    for gene in genes:
        if gene == "":
            pass
        else:
            if coordinates[i].split(" ")[0] == "None" or coordinates[i].split(" ")[0] == "N":
                position="None"
                contig="None"
            else:
                position=coordinates[i].split(" ")[0][1:]+"->"+coordinates[i].split(" ")[2][:-1]
                contig=coordinates[i].split(" ")[3].rstrip()
            salida.write(gene+"\t"+position+"\t"+contig+"\t"+pathway+"\t"+pathway_name+"\t"+microClass+"\t"+macroClass+"\n")
            i=i+1

