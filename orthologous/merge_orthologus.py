trato_list_file= open("ortologos.csv", "r")
trato_list = trato_list_file.readlines()
trato_list_file.close()
tatro_ids_jgi_dic={}
for line in trato_list:
    line_split=line.split(",")
    tatro_id=line_split[0]+"-T1"
    JGI_id=line_split[1]
    tatro_ids_jgi_dic[tatro_id]=JGI_id

def get_tatro_with_id(file_name):
    file = open(file_name, "r")
    lines = file.readlines()
    file.close()
    tatro_with_id={}
    for line in lines:
        line_split=line.split("\t")
        tatro_id=line_split[1]
        specie_id=line_split[2].rstrip()

        if "," in tatro_id:
            tatro_ids=tatro_id.split(",")
            for tatro_id in tatro_ids:
                tatro_id=tatro_id.replace(" ", "")
                tatro_with_id[tatro_id]=specie_id
        else:
            tatro_with_id[tatro_id]=specie_id
    
    return tatro_with_id

asper = get_tatro_with_id("Trichoderma_atroviride_IMI206040__v__Aspergillus_nidulans_ASM1142.tsv")
neuro = get_tatro_with_id("Trichoderma_atroviride_IMI206040__v__Neurospora_crassa_NC12.tsv")
T_harz = get_tatro_with_id("Trichoderma_atroviride_IMI206040__v__Trichoderma_harzianum_Tr1.tsv")
T_reseei = get_tatro_with_id("Trichoderma_atroviride_IMI206040__v__Trichoderma_reesei_QM6a.tsv")
T_virens = get_tatro_with_id("Trichoderma_atroviride_IMI206040__v__Trichoderma_virens_Gv29-8.tsv")

all_tatros_ids = list(tatro_ids_jgi_dic) + list(asper) + list(neuro) + list(T_harz) + list(T_reseei) + list(T_virens)
all_tatros_ids = (list(set(all_tatros_ids)))

output_file=open("ortologos_especies.csv", "w")

output_file.write("Atroviride,Tatro_JGI,Neurospora,Aspergillus,Harzianum,Reesei,Virens")
output_file.write("\n")

def to_output(dic, tatro_id):
    return dic.get(tatro_id,"-").replace(",", "")

for tatro_id in all_tatros_ids:
    output_string = ",".join([
     tatro_id,
     to_output(tatro_ids_jgi_dic, tatro_id), 
     to_output(asper, tatro_id), 
     to_output(neuro, tatro_id),
     to_output(T_harz, tatro_id), 
     to_output(T_reseei, tatro_id),
     to_output(T_virens, tatro_id)])

    output_file.write(output_string)
    output_file.write("\n")

output_file.close()