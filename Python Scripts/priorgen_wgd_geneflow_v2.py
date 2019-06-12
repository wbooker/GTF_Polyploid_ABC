#!/usr/bin/python
# -*- coding: utf-8 -*-
from numpy.random import uniform
from numpy.random import exponential
from numpy.random import beta
from numpy.random import binomial
from random import shuffle
from math import log
import sys
n1, n2, nA, tau, M, exp, alpha, shape1, shape2 = [], [], [], [], [], [], [], [], []
help="\n\t\033[32mExternal required library: numpy \033[1;m(sudo apt-get install python-numpy)\n\
\tpriorgen.py generates prior distributions for multiple multilocus simulations under 27 different models of speciation with/withou interploidal gene-flow. The output can be used from the stdout by ms (Hudson 2002), msnsam (Ross-Ibarra 2008) and msms (Ewing and Hermisson 2010) using the 'tbs' feature.\n\
\tIt requires one input file containing six lines: \n\
\t\tL1=description line, non read by priorgen.py\n\
\t\tL2=a vector with the lengths (L) for each of the surveyed locus\n\
\t\tL3=a vector with the number of sampled individuals (nspA) for each locus, for the first population\n\
\t\tL4=a vector with the number of sampled individuals (nspB) for each locus, for the second population\n\
\t\tL5=a vector with the populational mutation rates theta(i)=4.N.µ.L(i) for each locus 'i'\n\
\t\tL6=a vector with the populational recombination rates rho(i)=4.N.r.L(i) for each locus 'i'\n\n\
\tValues print in the stdout are used by ms-like coalescent simulators, values written in a file are the multilocus parameters useful for an ABC analysis\n\n\
\t\tparameters: name of the output file name. Ex \033[1;35mparameters=listOfParameters.txt\033[1;m\n\
\t\tn1: prior for N1 (the effective population size of the first population). Ex \033[1;35mn1=0 n1=10\033[1;m\n\
\t\tn2: prior for N2 (the effective population size of the second population). Ex\033[1;35m n2=0 n2=10\033[1;m\n\
\t\tnA: prior for NA (the effective population size of the ancetral population). Ex\033[1;35m nA=0 nA=10\033[1;m\n\
\t\ttau: prior for Tsplit (the time of speciation). Ex\033[1;35m tau=0 tau=3\033[1;m\n\
\t\talpha: prior for the proportion of tetrasomic genome, between 0 and 1, but can be fixed. Ex\033[1;35m alpha=0 alpha=1\033[1;m\n\
\t\tmigModel: model for interploidal gene-flow. \033[1;35mnoMig\033[1;m (no-migration), \033[1;35mmigA\033[1;m (migration between diploid and one sub-genome) or \033[1;35mmigAB\033[1;m (gene flow between diploid and the two sub-genomes)\n\
\t\tM: prior for migration rate 4.N.m. Ex\033[1;35m M=0 M=4\033[1;m\n\
\t\tshape1: prior for the first shape parameter of the Beta distribution. Ex\033[1;35m shape1=0 shape1=10\033[1;m\n\
\t\tshape2: prior for the second shape parameter of the Beta distribution. Ex\033[1;35m shape2=0 shape2=50\033[1;m\n\
\t\tmodel: \033[1;35m=auto\033[1;m (Autopolyploidization), \033[1;35m=allo\033[1;m (Allopolyploidization), or \033[1;35m=polyP\033[1;m (Polyploid speciation)\n\
\t\tMVariation: \033[1;35m=homo\033[1;m (shared values of M throughout genome) or \033[1;35m=hetero\033[1;m (variation of M throughout genome)\n\
\t\tnreps: number of multilocus simulations. Ex \033[1;35mnreps=1000\033[1;m\n\
\t\tsymMig: \033[1;35m=sym\033[1;m (diplo->tetra = tetra->diplo) or \033[1;35m=asym\033[1;m (diplo->tetra != tetra->diplo)\n\n\
\t\toneWay: Is migration only from diploid to tetraploid? Ex \033[1;35oneWay=T\033[1;m\n\
\t\tTauDist: exponential distribution or uniform for tau? If exp, uses 2nd Tau for mean, and 1st Tau for min(for Twgd only, recommend setting to 0) Ex \033[1;35TauDist=E, TauDist = U\033[1;m\n\
\tEx:\n\
\t\033[1;32m./priorgen_v2.py bpfile=bpfile_test.txt n1=0 n1=1 n2=1 n2=2 nA=2 nA=3 tau=0 tau=4 alpha=0 alpha=1 nreps=10 M=0 M=10 shape1=0 shape1=10 shape2=0 shape2=100 model=allo migModel=migAB symModel=asym MVariation=hetero oneWay=T parameters=output.txt\033[1;m\n\n\
\tmsnsam tbs 40 -t tbs -r tbs tbs -I 3 tbs tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 2 3 tbs -m 3 2 tbs -n 1 tbs -n 2 tbs -n 3 tbs -g 2 tbs -g 3 tbs -ej tbs 3 1 -en tbs 2 tbs -em tbs 1 2 0 -em tbs 2 1 0 -ej tbs 2 1 -eN tbs tbs\t#for 'allo'\n\
\tmsnsam tbs 40 -t tbs -r tbs tbs -I 3 tbs tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 2 3 tbs -m 3 2 tbs -n 1 tbs -n 2 tbs -n 3 tbs -g 2 tbs -g 3 tbs -ej tbs 3 2 -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs\t#for 'auto'\n\
\tmsnsam tbs 40 -t tbs -r tbs tbs -I 3 tbs tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 2 3 tbs -m 3 2 tbs -n 1 tbs -n 2 tbs -n 3 tbs -g 2 tbs -g 3 tbs -ej tbs 3 2 -ej tbs 2 1 -eN tbs tbs\t#for 'polyP'\n\n\
\tcamille.roux.1983@gmail.com\n\
\t04/12/2014\n"
for i in sys.argv:
	if("help" in i):
		print(help)
		sys.exit(0)

if len(sys.argv)<=1:
	print(help)
	sys.exit(0)

for i in sys.argv:
	if "=" in i:
		i=i.split("=")
		if(i[0]=="bpfile"):
			bpfile=i[1]
		if(i[0]=="parameters"):
			outputParameters=i[1]
		if(i[0]=="n1"):
			n1.append(float(i[1]))
		if(i[0]=="n2"):
			n2.append(float(i[1]))
		if(i[0]=="nA"):
			nA.append(float(i[1]))
		if(i[0]=="tau"):
			tau.append(float(i[1]))
		if(i[0]=="M"):	#taux de migration
			M.append(float(i[1]))
		if(i[0]=="model"):	#auto; allo; polyP
			model=i[1]
		if(i[0]=="nreps"):
			nreps=int(i[1])
		if(i[0]=="migModel"):	#noMig; MigA; MigAB
			migModel=i[1]
		if(i[0]=="symModel"):	#sym; asym
			symModel=i[1]
		if(i[0]=="MVariation"):	#homo; hetero
			MVariation=i[1]
		if(i[0]=="alpha"):	#alpha de la recombinaison
			alpha.append(float(i[1]))
                if(i[0]=="shape1"):
                        shape1.append(float(i[1]))
                if(i[0]=="shape2"):
                        shape2.append(float(i[1]))
		if(i[0]=="oneWay"):
			oneWay=i[1]
		if(i[0]=="TauDist"):
			TauDist = i[1]

def binomBeta(nlocus, shape1, shape2, scalar):
        neutre=[0]
        hetero=[1]
        nNeutre=int(uniform(0, nlocus))
        nHetero=nlocus-nNeutre
        status=nNeutre*neutre+nHetero*hetero
        shuffle(status)
        values=[]
        for i in status:
                if i==0:
                        values.append(scalar)
                if i==1:
                        values.append(scalar*beta(shape1, shape2))
        res={}
        res["values"]=values
        res["nNeutre"]=nNeutre
        return(res)

def inheritance(nlocus, alpha):
        diso=[0]
        tetra=[1]
        nTetra=int(binomial(nlocus, alpha))
        nDiso=nlocus-nTetra
        status=nTetra*tetra+nDiso*diso
        shuffle(status)
        values=[]
        for i in status:
                if i==0:
                        values.append(0)
                if i==1:
                        values.append(100)
        res={}
        res["values"]=values
        res["nTetra"]=nTetra
        return(res)

infile=open(bpfile, "r")
tmp=infile.readline()	#skip the header
L=infile.readline().strip().replace(" ", "\t")
L=L.split("\t")
nspA=infile.readline().strip().replace(" ", "\t")
nspA=nspA.split("\t")
nspB=infile.readline().strip().replace(" ", "\t")
nspB=nspB.split("\t")
theta=infile.readline().strip().replace(" ", "\t")
theta=theta.split("\t")
rho=infile.readline().strip().replace(" ", "\t")
rho=rho.split("\t")

nlocus=int(len(L))

res="N1\tN2\tN2A\tNA\tTsplit\tTwgd\tnTetrasomic\tM1A\tM1B\tMA1\tMB1\tshape1Mig1A\tshape1MigA1\tshape1Mig1B\tshape1MigB1\tshape2Mig1A\tshape2MigA1\tshape2Mig1B\tshape2MigB1\talphaM1A\talphaMA1\talphaM1B\talphaMB1\n"

for i in range(nreps):
	n1prior=uniform(n1[0], n1[1])
	n2prior=uniform(n2[0], n2[1])
	n2aprior=uniform(n2[0], n2[1])
	nAprior=uniform(nA[0], nA[1])
	if(TauDist=="U"):
		Tsplit=uniform(tau[0], tau[1])
	if(TauDist=="E"):
		Tsplit=exponential(tau[1])
	if(model=="allo" or model=="auto"):
		Twgd=uniform(min(tau), Tsplit)
	if(model=="polyP"):
		Twgd=Tsplit
	growth=-log(2/(100000*n2prior))*1/Twgd
	alphaprior=uniform(alpha[0], alpha[1])  #le Alpha de la recombinaison: Si=0, alors disomy. Si=1 alors tetrasomy
	recomb=inheritance(nlocus, alphaprior)
	recombGenomic=recomb["values"]
	nTetra=recomb["nTetra"]
	#dealing now with interploidal gene-flow
	# A->1 = B->1
	if (oneWay=="T"):
		M1Aprior=uniform(0, 0)
	else:
		M1Aprior=uniform(M[0], M[1])
	M1Bprior=M1Aprior	#uniform(M[0], M[1])
	shape1Mig1Aprior=uniform(shape1[0], shape1[1])
	shape1Mig1Bprior=shape1Mig1Aprior	#uniform(shape1[0], shape1[1])
	shape2Mig1Aprior=uniform(shape2[0], shape2[1])
	shape2Mig1Bprior=shape2Mig1Aprior	#uniform(shape2[0], shape2[1])
	# 1->1 = 1->B
	MA1prior=uniform(M[0], M[1])
	MB1prior=MA1prior	#uniform(M[0], M[1])
	shape1MigA1prior=uniform(shape1[0], shape1[1])
	shape1MigB1prior=shape1MigA1prior	#uniform(shape1[0], shape1[1])
	shape2MigA1prior=uniform(shape2[0], shape2[1])
	shape2MigB1prior=shape2MigA1prior	#uniform(shape2[0], shape2[1])
	if(migModel=="noMig"):
		#set all migration parameters to 0
		M1ApriorGenomic={}
		M1ApriorGenomic["values"]=[0]*nlocus
		MA1priorGenomic={}
		MA1priorGenomic["values"]=[0]*nlocus
		M1BpriorGenomic={}
		M1BpriorGenomic["values"]=[0]*nlocus
		MB1priorGenomic={}
		MB1priorGenomic["values"]=[0]*nlocus
		MA1prior=0
		shape1MigA1prior=-9
		shape2MigA1prior=-9
		alphaMA1=-9
		M1Aprior=0
		shape1Mig1Aprior=-9
		shape2Mig1Aprior=-9
		alphaM1A=-9
		MB1prior=0
		shape1MigB1prior=-9
		shape2MigB1prior=-9
		alphaMB1=-9
		M1Bprior=0
		shape1Mig1Bprior=-9
		shape2Mig1Bprior=-9
		alphaM1B=-9
	if(migModel=="migA"):
		# B->1 = 0
		M1BpriorGenomic={}
		M1BpriorGenomic["values"]=[0]*nlocus
		M1Bprior=0
		shape1Mig1Bprior=-9
		shape2Mig1Bprior=-9
		alphaM1B=-9
		# 1->B = 0
		MB1priorGenomic={}
		MB1priorGenomic["values"]=[0]*nlocus
		MB1prior=0
		shape1MigB1prior=-9
		shape2MigB1prior=-9
		alphaMB1=-9
		#set migration parameters between pop-1 and pop-A to random values
		# A -> 1
		#M1Aprior=uniform(M[0], M[1])
		# 1 -> A
		#MA1prior=uniform(M[0], M[1])
		if(MVariation=="hetero"):
			# A -> 1
			M1ApriorGenomic=binomBeta(nlocus=nlocus, shape1=shape1Mig1Aprior, shape2=shape2Mig1Aprior, scalar=M1Aprior)
			alphaM1A=M1ApriorGenomic["nNeutre"]     #alpha de la migration
			if(symModel=="sym"):
				# 1->A = A->1
				MA1priorGenomic={}
				MA1priorGenomic["values"]=M1ApriorGenomic["values"]
				MA1prior=M1Aprior
				shape1MigA1prior=shape1Mig1Aprior
				shape2MigA1prior=shape2Mig1Aprior
				alphaMA1=alphaM1A
			if(symModel=="asym"):
				# 1->A != A->1
				MA1priorGenomic=binomBeta(nlocus=nlocus, shape1=shape1MigA1prior, shape2=shape2MigA1prior, scalar=MA1prior)
				alphaMA1=MA1priorGenomic["nNeutre"]
		if(MVariation=="homo"):
			# A -> 1
			M1ApriorGenomic={}
			M1ApriorGenomic["values"]=[M1Aprior]*nlocus
			alphaM1A=nlocus
			shape1Mig1Aprior=-9
			shape2Mig1Aprior=-9
			if(symModel=="sym"):
				# A->1 = 1->A
				MA1priorGenomic={}
				MA1priorGenomic["values"]=M1ApriorGenomic["values"]
				MA1prior=M1Aprior
				shape1MigA1prior=shape1Mig1Aprior
				shape2MigA1prior=shape2Mig1Aprior
				alphaMA1=nlocus
			if(symModel=="asym"):
				# A->1 != 1->A
				MA1priorGenomic={}
				MA1priorGenomic["values"]=[MA1prior]*nlocus
				alphaMA1=nlocus
	if(migModel=="migAB"):
		#M1Aprior=uniform(M[0], M[1])
		#MA1prior=uniform(M[0], M[1])
		#M1Bprior=M1Aprior	#migration A->1 = migration B->1
		#MB1prior=MA1prior	#migration 1->A = migration 1->B
		if(MVariation=="hetero"):
			# A -> 1
			M1ApriorGenomic=binomBeta(nlocus=nlocus, shape1=shape1Mig1Aprior, shape2=shape2Mig1Aprior, scalar=M1Aprior)
			alphaM1A=M1ApriorGenomic["nNeutre"]     #alpha de la migration
			if(symModel=="sym"):
				# 1 -> A
				MA1priorGenomic={}
				MA1priorGenomic["values"]=M1ApriorGenomic["values"]
				MA1prior=M1Aprior
				shape1MigA1prior=shape1Mig1Aprior
				shape2MigA1prior=shape2Mig1Aprior
				alphaMA1=alphaM1A
				# 1 -> B
				MB1priorGenomic={}
				MB1priorGenomic["values"]=M1ApriorGenomic["values"]
				MB1prior=M1Aprior
				shape1MigB1prior=shape1Mig1Aprior
				shape2MigB1prior=shape2Mig1Aprior
				alphaMB1=alphaM1A
				# B -> 1
				M1BpriorGenomic={}
				M1BpriorGenomic["values"]=M1ApriorGenomic["values"]
				M1Bprior=M1Aprior
				shape1Mig1Bprior=shape1Mig1Aprior
				shape2Mig1Bprior=shape2Mig1Aprior
				alphaM1B=alphaM1A
			if(symModel=="asym"):
				# B->1 = A->1
				M1BpriorGenomic=M1ApriorGenomic
				alphaM1B=M1ApriorGenomic["nNeutre"]
				shape1Mig1Bprior=shape1Mig1Aprior
				# 1->A != A->1
				MA1priorGenomic=binomBeta(nlocus=nlocus, shape1=shape1MigA1prior, shape2=shape2MigA1prior, scalar=MA1prior)
				alphaMA1=MA1priorGenomic["nNeutre"]
				# 1->B = 1->A
				MB1priorGenomic=MA1priorGenomic
				alphaMB1=MA1priorGenomic["nNeutre"]
				shape1MigB1prior=shape1MigA1prior
		if(MVariation=="homo"):
			# A -> 1
			alphaMA1=nlocus
			alphaM1A=nlocus
			alphaMB1=nlocus
			alphaM1B=nlocus
			shape1MigA1prior=-9
			shape1Mig1Aprior=-9
			shape1MigB1prior=-9
			shape1Mig1Bprior=-9
			shape2MigA1prior=-9
			shape2Mig1Aprior=-9
			shape2MigB1prior=-9
			shape2Mig1Bprior=-9
			M1ApriorGenomic={}
			M1ApriorGenomic["values"]=[M1Aprior]*nlocus
			if(symModel=="sym"):
				# 1->A = A->1
				MA1priorGenomic={}
				MA1priorGenomic["values"]=M1ApriorGenomic["values"]
				MA1prior=M1Aprior
				# 1->B = A->1
				MB1priorGenomic={}
				MB1priorGenomic["values"]=M1ApriorGenomic["values"]
				MB1prior=M1Aprior
				# B->1 = A->1
				M1BpriorGenomic={}
				M1BpriorGenomic["values"]=M1ApriorGenomic["values"]
				M1Bprior=M1Aprior
			if(symModel=="asym"):
				# B->1 = A->1
				M1BpriorGenomic={}
				M1BpriorGenomic["values"]=[M1Aprior]*nlocus
				M1Bprior=M1Aprior
				# 1->A != A->1
				MA1priorGenomic={}
				MA1priorGenomic["values"]=[MA1prior]*nlocus
				# 1->B = 1->A
				MB1priorGenomic={}
				MB1priorGenomic["values"]=[MA1prior]*nlocus
				MB1prior=MA1prior
	res+="{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\n".format(n1prior, n2prior, n2aprior, nAprior, Tsplit, Twgd, nTetra, M1Aprior, M1Bprior, MA1prior, MB1prior, shape1Mig1Aprior, shape1MigA1prior, shape1Mig1Bprior, shape1MigB1prior, shape2Mig1Aprior, shape2MigA1prior, shape2Mig1Bprior, shape2MigB1prior, alphaM1A, alphaMA1, alphaM1B, alphaMB1)
	for loc in range(nlocus):
		cout=""
		if(model=="auto"):
			#msnsam int(nspA[loc])+int(nspB[loc]) 85200000 -t float(theta[loc]) -r float(rho[loc]) int(L[loc]) -I 3 int(nspA[loc]) int(nspB[loc]) int(nspB[loc]) 0 -m 1 2 M1BpriorGenomic["values"][loc] -m 2 1 MB1priorGenomic["values"][loc] -m 1 3 M1ApriorGenomic["values"][loc] -m 3 1 MA1priorGenomic["values"][loc] -m 2 3 recombGenomic[loc] -m 3 2 recombGenomic[loc] -n 1 n1prior -n 2 n2prior -n 3 n2prior -g 2 growth -g 3 growth -ej Twgd 3 2 -en Twgd 2 n2aprior -ej Tsplit 2 1 -eN Tsplit nAprior
			cout+="{0} {1:.5f} {2:.5f} {3} {4} {5} {6} {7:.5f} {8:.5f} {9:.5f} {10:.5f} {11:.5f} {12:.5f} {13:.5f} {14:.5f} {15:.5f} {16:.5f} {17:.5f} {18:.5f} {19:.5f} {20:.5f} {21:.5f} {22:.5f} {23:.5f}".format(int(nspA[loc])+int(nspB[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), int(nspB[loc]), M1BpriorGenomic["values"][loc], MB1priorGenomic["values"][loc], M1ApriorGenomic["values"][loc], MA1priorGenomic["values"][loc], recombGenomic[loc], recombGenomic[loc], n1prior, n2prior, n2prior, growth, growth, Twgd, Twgd, n2aprior, Tsplit, Tsplit, nAprior)
		if(model=="allo"):
			#msnsam int(nspA[loc])+int(nspB[loc]) 85200000 -t float(theta[loc]) -r float(rho[loc]) int(L[loc]) -I 3 int(nspA[loc]) int(nspB[loc]) int(nspB[loc]) 0 -m 1 2 M1BpriorGenomic["values"][loc] -m 2 1 MB1priorGenomic["values"][loc] -m 1 3 M1ApriorGenomic["values"][loc] -m 3 1 MA1priorGenomic["values"][loc] -m 2 3 recombGenomic[loc] -m 3 2 recombGenomic[loc] -n 1 n1prior -n 2 n2prior -n 3 n2prior -g 2 growth -g 3 growth -ej Twgd 3 1 -en Twgd 2 n2aprior -em Twgd 1 2 -em Twgd 2 1 -ej Tsplit 2 1 -eN Tsplit nAprior
                        cout+="{0} {1:.5f} {2:.5f} {3} {4} {5} {6} {7:.5f} {8:.5f} {9:.5f} {10:.5f} {11:.5f} {12:.5f} {13:.5f} {14:.5f} {15:.5f} {16:.5f} {17:.5f} {18:.5f} {19:.5f} {20:.5f} {21:.5f} {22:.5f} {23:.5f} {24:.5f} {25:.5f}".format(int(nspA[loc])+int(nspB[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), int(nspB[loc]), M1BpriorGenomic["values"][loc], MB1priorGenomic["values"][loc], M1ApriorGenomic["values"][loc], MA1priorGenomic["values"][loc], recombGenomic[loc], recombGenomic[loc], n1prior, n2prior, n2prior, growth, growth, Twgd, Twgd, n2aprior, Twgd, Twgd, Tsplit, Tsplit, nAprior)
		if(model=="polyP"):
			#msnsam int(nspA[loc])+int(nspB[loc]) 85200000 -t float(theta[loc]) -r float(rho[loc]) int(L[loc]) -I 3 int(nspA[loc]) int(nspB[loc]) int(nspB[loc]) 0 -m 1 2 M1BpriorGenomic["values"][loc] -m 2 1 MB1priorGenomic["values"][loc] -m 1 3 M1ApriorGenomic["values"][loc] -m 3 1 MA1priorGenomic["values"][loc] -m 2 3 recombGenomic[loc] -m 3 2 recombGenomic[loc] -n 1 n1prior -n 2 n2prior -n 3 n2prior -g 2 growth -g 3 growth -ej Twgd 3 2 -ej Tsplit 2 1 -eN Tsplit nAprior
                        cout+="{0} {1:.5f} {2:.5f} {3} {4} {5} {6} {7:.5f} {8:.5f} {9:.5f} {10:.5f} {11:.5f} {12:.5f} {13:.5f} {14:.5f} {15:.5f} {16:.5f} {17:.5f} {18:.5f} {19:.5f} {20:.5f} {21:.5f}".format(int(nspA[loc])+int(nspB[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), int(nspB[loc]), M1BpriorGenomic["values"][loc], MB1priorGenomic["values"][loc], M1ApriorGenomic["values"][loc], MA1priorGenomic["values"][loc], recombGenomic[loc], recombGenomic[loc], n1prior, n2prior, n2prior, growth, growth, Twgd, Tsplit, Tsplit, nAprior)
		print(cout)
outputfile=open(outputParameters, "w")
outputfile.write(res)
outputfile.close()
