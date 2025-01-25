from filesParse import *

blast2go_omicsbox= b2GoTable("blast2go_omicsbox_table.txt", 2, 9)
blast2go_gff_export = b2GoGffExport("blast2go_gff_export.gff")
blast2go_annot= b2GoAnnot("blast2go.annot")
eggnog_go_slim = b2GoTable("eggnog_trichoderma_atroviride_imi206040_v3_proteins_fa_eggnog_go_slim.txt", 0, 2)
eggnog=b2GoTable("eggnog_Trichoderma_atroviride_IMI206040_V3.proteins.fa_eggnog.txt", 1, 9)
iprscan=b2GOIprcscn("iprscan_Tatro_V3_annotations.txt")


def goTermsList(dictionary):
    GO_list=[]
    for gen in dictionary:
        GO_list = GO_list + dictionary[gen]
    return list(set(GO_list))


blast2go_omicsbox_gen=list(blast2go_omicsbox.keys())
blast2go_omicsbox_goTerms = goTermsList(blast2go_omicsbox)

blast2go_gff_export_gen=list(blast2go_gff_export.keys())
blast2go_gff_export_goTerms = goTermsList(blast2go_gff_export)

blast2go_annot_gen=list(blast2go_annot.keys())
blast2go_annot_goTerms = goTermsList(blast2go_annot)

eggnog_go_slim_gen=list(eggnog_go_slim.keys())
eggnog_go_slim_goTerms = goTermsList(eggnog_go_slim)

eggnog_gen=list(eggnog.keys())
eggnog_goTerms = goTermsList(eggnog)

iprscan_gen=list(iprscan.keys())
iprscan_goTerms = goTermsList(iprscan)


from supervenn import supervenn
import matplotlib.pyplot as plt

sets = [set(blast2go_omicsbox_gen), 
        set(blast2go_gff_export_gen),
        set(blast2go_annot_gen),
        set(eggnog_go_slim_gen),
        set(eggnog_gen),
        set(iprscan_gen),
        ]

labels =['blast2go_omicsbox_gen',
         'blast2go_gff_export_gen',
         'blast2go_annot_gen',
         'eggnog_go_slim_gen',
         'eggnog_gen',
         'iprscan_gen',
         ]

# plt.figure(figsize=(20, 10))
# supervenn(sets, labels, rotate_col_annotations=True,
#           col_annotations_area_height=1.2, sets_ordering='minimize gaps',
#           min_width_for_annotation=180)
# plt.show()

# supervenn(sets, labels, side_plots=False, widths_minmax_ratio=0.05)
# plt.show()

sets = [set(blast2go_omicsbox_goTerms), 
        set(blast2go_gff_export_goTerms),
        set(blast2go_annot_goTerms),
        set(eggnog_go_slim_goTerms),
        set(eggnog_goTerms),
        set(iprscan_goTerms)
        ]

labels =['blast2go_omicsbox_goTerms',
         'blast2go_gff_export_goTerms',
         'blast2go_annot_goTerms',
         'eggnog_go_slim_goTerms',
         'eggnog_goTerms',
         'iprscan_goTerms'
         ]

plt.figure(figsize=(20, 10))
supervenn(sets, labels, rotate_col_annotations=True,
          col_annotations_area_height=1.2, sets_ordering='minimize gaps',
          min_width_for_annotation=180)
plt.show()

supervenn(sets, labels, side_plots=False, widths_minmax_ratio=0.05)
plt.show()
