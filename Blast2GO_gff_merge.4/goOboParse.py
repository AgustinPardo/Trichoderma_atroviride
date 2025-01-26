from goatools.obo_parser import GODag

obodag=GODag("inputs/go.obo", optional_attrs={'consider', 'replaced_by'}, load_obsolete=True, prt=None) 

def goOboNames(go_id):
	go_name = obodag[go_id].name
	go_namespace = obodag[go_id].namespace
	return [[go_id],[go_name],[go_namespace]]

def goOboDeprecated(go_id):
	consider=obodag[go_id].consider
	if len(consider) != 0:
		consider_list=[]
		for term in consider:
			if term[:2] != "CL":
				consider_list=consider_list + [term]
			else:
				consider_list=consider_list + [(go_id,term)]
				
	# Include the deprecated?
	# consider_list=consider_list+[go_id]
	else:
		# Include or not the ones without consider option
		# consider_list=[]
		consider_list=[go_id]
	return consider_list

def goObo(lista_go_id):
	go_list_depr=[]
	go_list=[]
	for go_id in lista_go_id:
		if obodag[go_id].is_obsolete:
			go_list_depr=goOboDeprecated(go_id)+go_list_depr
		else:
			go_list=[go_id]+go_list

	go_list_ND = go_list+go_list_depr

	go_list_id=[]
	go_list_name=[]
	go_list_namespace=[]
	
	for go in go_list_ND:
		if type(go) is tuple:
			go_list_id = goOboNames(go[0])[0]+go_list_id
			go_list_name = goOboNames(go[0])[1]+go_list_name
			go_list_namespace = goOboNames(go[0])[2]+go_list_namespace
		else:
			go_list_id = goOboNames(go)[0]+go_list_id
			go_list_name = goOboNames(go)[1]+go_list_name
			go_list_namespace = goOboNames(go)[2]+go_list_namespace
	return go_list_id, go_list_name, go_list_namespace