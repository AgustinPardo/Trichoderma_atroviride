import Bio.SeqIO as bpsio 
import Bio.SeqIO as bpio
import sys

allcontigs = list(bpio.parse("NCBI_GCF_000171015.1_TRIAT_v2.0_genomic.gb", "gb"))
NdeCont=1               
for contig in range(len(allcontigs)):    
	res = allcontigs[contig]  # tomo un contig.
   
	# Contadores para caracterizar los genebanks.
	i = 2
	ORFs=0
	rRNAs=0
	tRNAs=0
	hyp=0
	hypconser=0
	putative=0
	uncha=0
	conEc=0
	notes=0
	while i < (len(res.features)):  # Recorro features: source, gen, CDS, rRNA, tRNA.
		if res.features[i].type == "tRNA":
			tRNAs+=1
		if res.features[i].type == "rRNA":
			rRNAs+=1
		if res.features[i].type == "CDS":
			
			if "product" in res.features[i].qualifiers:
				producto= res.features[i].qualifiers["product"][0]
				productoSeparado= producto.split(" ")
				productoUnido="-".join(productoSeparado)
				B=""
			   
				if "hypothetical protein" in producto and "conserved" not in producto:
					hyp+=1
		  
				if "conserved hypothetical protein" in producto:
					hypconser+=1
				
				if "putative membrane protein" in producto:
					putative+=1                    
				if  "putative signal peptide protein"  in producto:
					putative+=1
				if  "putative lipoprotein" in producto:
					putative+=1
							  
				if "ncharacteri" in producto:
					uncha+=1
			else:
				if "note" in res.features[i].qualifiers:
					nota=res.features[i].qualifiers["note"][0]
					if "hypothetical protein" in nota and "conserved" not in nota:
						hyp+=1
					if "conserved hypothetical protein" in nota:
						hypconser+=1
					if "ncharacteri" in nota:
						uncha+=1
					notes+=1 
					B=""   
				  
			if "EC_number" in res.features[i].qualifiers:
				listaEC=res.features[i].qualifiers["EC_number"]
				for x in range(len(listaEC)):
					if x==0:
						B=B+listaEC[x]
					else:
						B=B+","+listaEC[x]

				conEc+=1
			else:
				pass
			i = i + 1
			ORFs+=1
			
		else:
			i = i + 1
	NdeCont+=1            
	print("Contig "+str(contig))
	print("ORFS "+str(ORFs))
	print("tRNAs"+str(tRNAs))
	print("RNAs"+str(rRNAs))   
	print("hypothetical protein "+str(hyp))
	print("hypothetical conserved "+str(hypconser))       
	print("putative "+str(putative))
	print("uncharacterized "+str(uncha))
	print("ECs: "+str(conEc))
	print("note "+str(notes)) 
