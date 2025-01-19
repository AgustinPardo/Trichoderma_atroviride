def goCounts(dic_input):
	"""
	Method that describe the relation of gens and GO terms
    Args:
        dic_input: Dictiornary of gen(key) and GO terms as list(value)
    Returns:
        Print some basics stats of the relation of gene a Go terms"
	"""
	gen_conGO=0
	gen_sinGO=0
	gen_conmasdeunGO=0
	GO_totales=[]

	for gene in dic_input:
		go_terms_list = list(set(dic_input[gene]))
		GO_totales = GO_totales + go_terms_list
		if len(go_terms_list) == 0:
			gen_sinGO=1+gen_sinGO
		else:
			if len(go_terms_list)==1:
				gen_conGO=1+gen_conGO
			else:
				gen_conmasdeunGO=1+gen_conmasdeunGO	

	print("	Genes with 1 GO", gen_conGO)
	print("	Genes with more than 1 GO", gen_conmasdeunGO)
	print("	Genes without GO", gen_sinGO)
	print("	Total genes", len(dic_input))
	print("	Total GO terms", len(list(set(GO_totales))))

def outputReview(func):
  """
  Decorator that apply goCounts to result of function decorated.
  Args:
    func: The function to be decorated.
  Returns:
    A function that wrapped the result of the decorated function with goCounts.
  """
  def wrapper(*args, **kwargs):
    result = func(*args, **kwargs)
    print(*args)
    goCounts(result)
    return result
  return wrapper