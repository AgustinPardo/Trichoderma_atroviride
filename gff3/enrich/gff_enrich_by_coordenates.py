from BCBio import GFF

main_file = "Trichoderma_atroviride_chroms.gff3"
lncRNA_enrich_source_file = "Trichoderma_atroviride_lncRNA_chroms.gff3"
sRNA_enrich_source_file = "Trichoderma_atroviride_sRNA_chroms.gff3"
out_put_file = "gff_enriched.gff3"

main_handle = open(main_file)
out_handle = open(out_put_file, "w")

def parse_dic_to_str(dictionary):
    res=[]
    for key in dictionary:
        res.append(key+":"+dictionary[key][0])
    return "|".join(res)

def enrich_feature(source_file, enrich_type, main_feature, main_rec_id):
        enrich_source_handle = open(source_file)
        main_loc=main_feature.location   
        for source_rec in GFF.parse(enrich_source_handle):
            source_rec_id = source_rec.id
            if main_rec_id==source_rec_id:
                for source_feature in source_rec.features:
                    source_loc=source_feature.location
                    if (main_loc.start<=source_loc.start<=main_loc.end) or (main_loc.start<=source_loc.end<=main_loc.end):
                        # print(
                        #     enrich_type,
                        #     main_rec_id,
                        #     main_feature.id,
                        #     main_loc,
                        #     "->",
                        #     source_rec_id,
                        #     source_feature.id,
                        #     source_loc
                        #     )
                        main_feature.qualifiers[enrich_type]=parse_dic_to_str(source_feature.qualifiers)
        enrich_source_handle.close()

main_gff = GFF.parse(main_handle)
for main_rec in main_gff:
    main_rec_id=main_rec.id
    for main_feature in main_rec.features:
        if main_feature.type=="gene":          
            enrich_feature(lncRNA_enrich_source_file, "lncRNA", main_feature, main_rec_id)
            enrich_feature(sRNA_enrich_source_file, "sRNA", main_feature, main_rec_id)
GFF.write(main_gff, out_handle)

main_handle.close()
out_handle.close()
