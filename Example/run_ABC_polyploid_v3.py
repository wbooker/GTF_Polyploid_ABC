#!/usr/bin/python
## ./run_ABC_polyploid.py [poly mode] [inheritance mode] [migration model] [T/F one way migration] [bpfile] [num sims per rep] [begin folder] [end folder] [T/F keep the jsfs]
## model = [diso1, diso2_auto, diso2_allo, tetra1, tetra2_auto, tetra2_allo, hetero1, hetero2_auto, hetero2_allo]
listModels = ["auto", "allo", "polyP"]

import sys
import os
import shutil
import time
import datetime

help = "\n\trun_ABC_polyploid.py is a python script usefull to run some random simulations of models\n\tof speciation between a diploid and a tetraploid\n"
help += "\tIt takes 5 arguments:\n\t\t1) a demographic model (auto, allo);\n\t\t2) a model of inheritance (disomic, heterosomic, tetrasomic);\n\t\t3) a migration model; \n\t\t4) the name of the bpfile; \n\t\t5) the number of replicates (an integer);\n\t\t5) the number of replicate folders (total sims = replicates*folders)\n"
help += "\n\tThe accepted polyploid models are:\n\t\t{0}\n".format("\n\t\t".join(listModels))
help += "\n\tThe accepted polyploid models are:\n\t\t noMig \n\t\t migA \n\t\t migAB"
help += "\n\tcommand line: ./run_ABC_polyploid.py [polyploid model] [inheritance] [migration model] [bpfile] [nreps] [nfolders]\n\n"

if len(sys.argv) != 10:
	print(help)
	print("\n\tThe expected number of arguments is 9\n\n")
	sys.exit()

model = sys.argv[1]
inheritance = sys.argv[2]
migmodel = sys.argv[3]
oneWayIn = sys.argv[4]
bpfile = sys.argv[5]
nreps = int(sys.argv[6])
start = int(sys.argv[7])
nfolders = int(sys.argv[8])
jsfsTF = sys.argv[9]




if model not in listModels:
	print(help)
	print("\n\tThe model {0} is not in the list of expected models:\n{1}\n".format(model, "\n\t".join(listModels)))
	sys.exit()

########### Set prior parameters

n1_1=0
n1_2=20
n2_1=10
n2_2=50
nA_1=0
nA_2=50
tau_1=0
tau_2=3
TauDist="E"
M_1=0
M_2=10
shape1_1=0
shape1_2=10
shape2_1=0
shape2_2=100
symModel="asym"
MVariation="hetero"

#####################################

if os.path.isfile(bpfile) == False:
	print(help)
	print("\n\t{0} is not found as a correct bpfile\n".format(bpfile))
	sys.exit()

if inheritance == 'disomic':
	alpha1=0
	alpha2=0

if inheritance == 'tetrasomic':
	alpha1=1
	alpha2=1

if inheritance == 'heterosomic':
	alpha1=0
	alpha2=1




# get the number of loci
infile = open(bpfile, 'r')
tmp = infile.readline()
tmp = infile.readline().strip().split('\t')
nLoci = len(tmp)
infile.close()

# get the number of simulations
#infile = open(argfile, 'r')
#for i in infile:
#	if 'nreps' in i:
#		i = i.strip().split('=')
#		if i[0].strip() == 'nreps':
#			nreps = int(i[1])
#infile.close()

nsim = nreps * nLoci
if(oneWayIn == "T"):
	modeldir = "{0}_{1}_{2}_{3}".format(model,inheritance,migmodel,"1W")
if(oneWayIn == "F"):
	modeldir = "{0}_{1}_{2}".format(model,inheritance,migmodel)
continue_pipeline = "T"

if os.path.exists(modeldir):
	deldir = raw_input("model directory already exists, delete it and its contents and continue? (y/n)")
	if deldir == "y":
		shutil.rmtree(modeldir)
		os.mkdir(modeldir)
	if deldir == "n":
		print ("not writing over directory, exiting \n")
		continue_pipeline = "F"
else:
	os.mkdir(modeldir)

if model == "auto":
	command = "./priorgen_wgd_geneflow_v2.py bpfile={0} n1={9} n1={10} n2={11} n2={12} nA={13} nA={14} tau={15} tau={16} TauDist={17} alpha={4} alpha={5} nreps={1} M={18} M={19} shape1={20} shape1={21} shape2={22} shape2={23} model={2} migModel={6} symModel={24} MVariation={25} oneWay={8} parameters={7}/output.txt | msnsam tbs {3} -t tbs -r tbs tbs -I 3 tbs tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 2 3 tbs -m 3 2 tbs -n 1 tbs -n 2 tbs -n 3 tbs -g 2 tbs -g 3 tbs -ej tbs 3 2 -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs | ./mscalc_wgd_v3.py {7} {0}".format(bpfile, nreps, model, nsim, alpha1, alpha2, migmodel, modeldir, oneWayIn, n1_1, n1_2, n2_1, n2_2, nA_1, nA_2, tau_1, tau_2, TauDist, M_1, M_2, shape1_1, shape1_2, shape2_1, shape2_2, symModel, MVariation)

if model == "allo":
	command = "./priorgen_wgd_geneflow_v2.py bpfile={0} n1={9} n1={10} n2={11} n2={12} nA={13} nA={14} tau={15} tau={16} TauDist={17} alpha={4} alpha={5} nreps={1} M={18} M={19} shape1={20} shape1={21} shape2={22} shape2={23} model={2} migModel={6} symModel={24} MVariation={25} oneWay={8} parameters={7}/output.txt | msnsam tbs {3} -t tbs -r tbs tbs -I 3 tbs tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 2 3 tbs -m 3 2 tbs -n 1 tbs -n 2 tbs -n 3 tbs -g 2 tbs -g 3 tbs -ej tbs 3 1 -en tbs 2 tbs -em tbs 1 2 0 -em tbs 2 1 0 -ej tbs 2 1 -eN tbs tbs | ./mscalc_wgd_v3.py {7} {0}".format(bpfile, nreps, model, nsim, alpha1, alpha2, migmodel, modeldir, oneWayIn, n1_1, n1_2, n2_1, n2_2, nA_1, nA_2, tau_1, tau_2, TauDist, M_1, M_2, shape1_1, shape1_2, shape2_1, shape2_2, symModel, MVariation)

if model == "polyP":
	command = "./priorgen_wgd_geneflow_v2.py bpfile={0} n1={9} n1={10} n2={11} n2={12} nA={13} nA={14} tau={15} tau={16} TauDist={17} alpha={4} alpha={5} nreps={1} M={18} M={19} shape1={20} shape1={21} shape2={22} shape2={23} model={2} migModel={6} symModel={24} MVariation={25} oneWay={8} parameters={7}/output.txt | msnsam tbs {3} -t tbs -r tbs tbs -I 3 tbs tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 2 3 tbs -m 3 2 tbs -n 1 tbs -n 2 tbs -n 3 tbs -g 2 tbs -g 3 tbs -ej tbs 3 2 -ej tbs 2 1 -eN tbs tbs	| ./mscalc_wgd_v3.py {7} {0}".format(bpfile, nreps, model, nsim, alpha1, alpha2, migmodel, modeldir, oneWayIn, n1_1, n1_2, n2_1, n2_2, nA_1, nA_2, tau_1, tau_2, TauDist, M_1, M_2, shape1_1, shape1_2, shape2_1, shape2_2, symModel, MVariation)

print("Total number of datasets (nreps*nfolders) = {0}".format(nreps*(nfolders-start+1)))


try:
	if continue_pipeline == "T":
		for i in range(start,nfolders+1):
			if os.path.exists("seedms"):
				os.remove("seedms")
			ts = time.time()
			timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
			workdir = "{0}/{0}_{1}".format(modeldir,i)
			print("Replicate: {0}/{1} Started: {2}".format(i,nfolders,timestamp))
			os.system(command)
			os.mkdir(workdir)
			shutil.move("{0}/output.txt".format(modeldir),workdir)
			shutil.move("{0}/ABCstat.txt".format(modeldir),workdir)
			if (jsfsTF=="T"):
				shutil.move("{0}/ABCjsfs.txt".format(modeldir),workdir)
			if (jsfsTF=="F"):
				os.remove("{0}/ABCjsfs.txt".format(modeldir))
			os.remove("seedms")
except NameError:
	print("I don't know exactly why, but it doesn't work, sorry...")
