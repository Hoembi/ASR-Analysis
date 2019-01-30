from pymol import cmd, stored, math

def loadBnoputty (mol,startaa=1,source="newBfactors.txt", visual="Y"):
	"""
	Replaces B-factors with a list of values contained in a plain txt file

	usage: loadBnoputty mol, [startaa, [source, [visual]]]

	mol = any object selection (within one single object though)
	startaa = number of first amino acid in 'new B-factors' file (default=1)
	source = name of the file containing new B-factor values (default=newBfactors.txt)
	visual = redraws and colors the structre and displays bar with min/max values (default=Y)
	"""
	obj=cmd.get_object_list(mol)[0]
	cmd.alter(mol,"b=-1.0")
	inFile = open(source, 'r')
	counter=int(startaa)
	bfacts=[]
	for line in inFile.readlines():
		bfact=float(line)
		bfacts.append(bfact)
		cmd.alter("%s and resi %s and n. CA"%(mol,counter), "b=%s"%bfact)
		counter=counter+1
	if visual=="Y":
		cmd.show_as("cartoon",mol)
		cmd.spectrum("b","rainbow", "%s and n. CA " %mol)
		cmd.ramp_new("count", obj, [min(bfacts), max(bfacts)], "rainbow")
		cmd.recolor()

cmd.extend("loadBnoputty_recolor", loadBnoputty);
