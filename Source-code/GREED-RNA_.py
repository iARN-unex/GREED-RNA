import sys
import time
import random
import RNA
from pathlib import Path
import argparse
from difflib import SequenceMatcher
from itertools import combinations

############################################################################################################################################################################################################################################
#       RNA STRUCTURE ANALYSIS
############################################################################################################################################################################################################################################
def getComponents_B_U(_struct, init=0):
    a=0
    b=0
    up_=[]
    U_ =[]
    B_ =[]
    for i in range(0,len(_struct)):
            if _struct[i] == '.':
                    up_.append(init+i)
            elif _struct[i] == '(':
                    a+=1
                    c=0
                    for j in range((i+1),len(_struct)):
                            if _struct[j] == '(':
                                    c+=1
                            elif _struct[j] == ')':
                                    c-=1
                                    if c == -1:
                                            U_.append(init+i)
                                            B_.append(init+j)
                                            break
            elif _struct[i] == ')':
                    b+=1
            else:
                    sys.exit('_structure error: Wrong character')
    if a != b:
            sys.exit('_structure error: Numbers of "(" and ")" do not match')

    U_ = U_ + up_

    return B_,U_;

def printComponents_B_U(B, U):
	print ("B: ", end='')
	for i in range(len(B)):
	        print ("(" + str(U[i]) + "," + str(B[i]) + ") ", end='')
	print ()
	print ("U: ", end='')
	for i in range(len(B),len(U)):
	        print (str(U[i]) + " ", end='')
	print ()

############################################################################################################################################################################################################################################
#       EVALUATION
############################################################################################################################################################################################################################################
def calculate_pGC(sequence, structure):
	B_,U_ = getComponents_B_U(structure)

	#basepairs distribution
	pairs=[]
	for i in range(0,len(B_)):
        	j = U_[i]
        	k = B_[i]
        	pairs.append(sequence[j]+sequence[k])

	CG=pairs.count('CG')/max(len(pairs),1)
	GC=pairs.count('GC')/max(len(pairs),1)
	pGC=CG+GC

	return pGC

def hammingDistance(str1, str2):
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

def calculateCompositionM2dRNAs(sequence, structure):
        B_,U_ = getComponents_B_U(structure)

        #basepairs distribution
        pairs=[]
        for i in range(0,len(B_)):
                j = U_[i]
                k = B_[i]
                pairs.append(sequence[j]+sequence[k])

        AU=pairs.count('AU')/max(len(pairs),1)
        UA=pairs.count('UA')/max(len(pairs),1)
        CG=pairs.count('CG')/max(len(pairs),1)
        GC=pairs.count('GC')/max(len(pairs),1)
        GU=pairs.count('GU')/max(len(pairs),1)
        UG=pairs.count('UG')/max(len(pairs),1)

        pAU=AU+UA
        pGC=CG+GC
        pUG=UG+GU
        paired = max(pAU,pGC,pUG)

        #unpaired bases distribution
        subseq=""
        for i in range(len(B_),len(U_)):
                j = U_[i]
                subseq=subseq+sequence[j]

        uA=subseq.count('A')/len(subseq)
        uU=subseq.count('U')/len(subseq)
        uG=subseq.count('G')/len(subseq)
        uC=subseq.count('C')/len(subseq)

        unpaired = max(uA,uU,uG,uC)

        #total bases distribution
        a=sequence.count('A')/len(sequence)
        u=sequence.count('U')/len(sequence)
        g=sequence.count('G')/len(sequence)
        c=sequence.count('C')/len(sequence)

        total = max(a,u,g,c)

        return round((total+unpaired+paired)*100,3), pGC

def evaluateM2DRNAS(s, t):
	# f1=partition_function
	# f2=ensemble_diversity
	# f3=composition
	# c1=similarity
	# RNA_SEQ
	# STRUCTURE

	fc = RNA.fold_compound(s)
	(mfe_struct, mfe) = fc.mfe()

	f_bpdistance = RNA.bp_distance(mfe_struct, t)

	fc.exp_params_rescale(mfe)
	(pp, f_pf) = fc.pf()

	f_ensemble = fc.mean_bp_distance()
	f_composition,pGC = calculateCompositionM2dRNAs (s, mfe_struct)

	return f_pf, f_ensemble, f_composition, f_bpdistance, s, mfe_struct, pGC

def evaluate(s, t):
	fc = RNA.fold_compound(s)
	(mfe_struct, mfe) = fc.mfe()
	f_bpdistance = RNA.bp_distance(mfe_struct, t)
	fc.exp_params_rescale(mfe)
	(pp, f_pf) = fc.pf()
	f_ensemble = fc.ensemble_defect(t)
	f_prob = fc.pr_structure(t)
	f_dist = hammingDistance(mfe_struct, t)
	f_pGC = calculate_pGC (s, mfe_struct)

	return  mfe_struct, round(f_bpdistance/len(t),4), (-1.0)*f_prob, round(f_pf,4), round(f_ensemble,4), round(f_pGC,4), round(f_dist/len(t),4)

############################################################################################################################################################################################################################################
#       GREEDY INITIALIZATION
############################################################################################################################################################################################################################################
def greedyInitialization(target, B, U):
    s = list('_'*len(target))
    triplet = random.choice(['GCA', 'CGA'])
    for i in range(len(B)):
      if s[B[i]] == '_' and s[U[i]] == '_':
        s[B[i]] = triplet[0]
        s[U[i]] = triplet[1]
    for i in range(len(B), len(U)):
      if s[U[i]] == '_':
        s[U[i]] = triplet[2]
    return "".join(s)

############################################################################################################################################################################################################################################
#       RANDOM MUTATION
############################################################################################################################################################################################################################################
def randomMutation (sequence, pm, B, U):
	s = list(sequence)
	for x in range (pm):
		pos = int(random.uniform(0.00, len(U)))
		if pos < len(B):
			i = U[pos]
			j = B[pos]
			pair = random.choice(['GC', 'CG', 'AU', 'UA', 'GU', 'UG'])
			s[i] = pair[0]
			s[j] = pair[1]
		else:
			i = U[pos]
			unpair = random.choice(['A' , 'C', 'G', 'U'])
			s[i] = unpair
	return "".join(s)

############################################################################################################################################################################################################################################
#       GREEDY MUTATION
############################################################################################################################################################################################################################################
def greedyMutation (sequence, pm, B, U):
	s = list(sequence)
	for x in range (pm):
		pos = int(random.uniform(0.00, len(U)))
		if pos < len(B):
			i = U[pos]
			j = B[pos]
			pair = random.choice(['GC', 'CG'])
			s[i] = pair[0]
			s[j] = pair[1]
		else:
			i = U[pos]
			unpair = random.choice(['A'])
			s[i] = unpair
	return "".join(s)

############################################################################################################################################################################################################################################
#       GREEDy Evolutionary Strategy
############################################################################################################################################################################################################################################
def greed (target, maxGC, maxN, pm, maxStagnation, stop):

	start = time.time()
	B,U = getComponents_B_U(target)

	totalRuntimeVienna = 0

	amount = max(round(len(target) * pm),1)
	sequences_found = [] # Exactly folds into target, maxGC is not taken into account at this point
	n_evals = 0
	stagnation = 0
	flagInit = False
	flagBPD = False
	flagGC = False
	randomMutationFlag = False
	minBPD = 999999
	minGC = 999999
	minN = 0

	parent_seq = ''
	parent_eval = []
	offspr_seq = ''
	offspr_eval = []

	iter = 0

	print()
	print()
	print (" (1) Let's minimize the base-pair distance between predicted and target structures")
	print ("    -> Min(Base-pair distance), min(Hamming distance), max(probability over ensemble), min(partition function), min(ensemble defect), min(%GC-content)")
	print()
	print()
	while time.time() - start < stop:
		iter += 1

		# GREEDY INITIALIZATION (first time and after re-start)
		if parent_seq == '':
			flagInit = True
			parent_seq = greedyInitialization(target, B, U)
			_initEvalTime = time.time()
			parent_eval = evaluate(parent_seq, target)
			totalRuntimeVienna += (time.time() - _initEvalTime)
			n_evals += 1

		# MUTATION GREEDY-RANDOM
		if randomMutationFlag == False:
			offspr_seq = greedyMutation(parent_seq, amount, B, U)
		else:
			offspr_seq = randomMutation(parent_seq, amount, B, U)

		# EVALUATE OFFSPRING
		if offspr_seq != parent_seq:
			_initEvalTime = time.time()
			offspr_eval = evaluate(offspr_seq, target)
			totalRuntimeVienna += (time.time() - _initEvalTime)
			n_evals += 1
		else:
			offspr_eval = parent_eval

		results = []
		results.append([parent_seq] + list(parent_eval))
		results.append([offspr_seq] + list(offspr_eval))

		# SELECTION DYNAMICALLY
		#  (0) seq
		#  (1) mfe_struct
		#  (2) f_bpdistance/len(t)
		#  (3) f_prob
		#  (4) f_p
		#  (5) f_ensemble
		#  (6) f_pGC
		#  (7) f_dist/len(t)
		#
		if flagBPD == True:
			if flagGC == False:
				results.sort(key=lambda x: (x[2], x[7], x[6], x[4], x[5], x[3]))
			else:
				results.sort(key=lambda x: (x[2], x[7], abs(x[6]-maxGC), x[4], x[5], x[3]))
		else:
			results.sort(key=lambda x: (x[2], x[7], x[3], x[5], x[4], x[6]))

		if parent_seq != results[0][0] or flagInit == True:

			if parent_eval[1] > results[0][2] or (parent_eval[1] == 0 and parent_eval[5] > results[0][6]) or flagInit == True:
				print ('\t ->' + str([parent_seq] + list(parent_eval)))

			flagInit = False
			stagnation = 0
			parent_seq = results[0][0]
			parent_eval = results[0][1:]

			if parent_eval[1] == 0:
				if flagBPD == False:
					flagBPD = True
					print ()
					print ()
					print ("    -> Found sequence! Elapsed time:                   " + str(time.time() - start_time) + " secs.")
					print ()
					print ()
					print (" (2) We found a sequence that folds exactly like the target structure. Let's minimize the GC-content")
					print ("    -> ORDER: Min(Base-pair distance), min(Hamming distance), min(GC-content), min(partition function), min(ensemble defect), max(probability over ensemble))")
					print ()
					print ()

				if parent_seq not in sequences_found:
					sequences_found.append(parent_seq)

				if parent_eval[5] <= maxGC:
					flagGC = True
					randomMutationFlag = True
					minN += 1
					if flagGC == False:
						print ()
						print ()
						print ("    -> Satisfied GC-constraint! Elapsed time (2):                   " + str(time.time() - start_time) + " secs.")
						print ()
						print ()
						print (" (3) We found 1 sequence that folds exactly like the target structure with a GC-content lower than " + str(maxGC) + ". Let's find " + str(maxN - minN) + " more structures")
						print ("    -> ORDER: Min(Base-pair distance), min(Hamming distance), min(abs(GC-content-maxGC)), min(partition function), min(ensemble defect), max(probability over ensemble))")
						print ()
						print ()

			if minBPD > parent_eval[1] or (parent_eval[1] == 0 and minGC > parent_eval[5]):
				minBPD = parent_eval[1]
				if minGC > parent_eval[5]:
					minGC = parent_eval[5]
		else:
			stagnation += 1
			if flagBPD == True and flagGC == True and minN >= maxN:
				print ()
				print ()
				print ("    -> Finished!! Elaapsed time:                   " + str(time.time() - start_time) + " secs.")
				break
			elif stagnation >= maxStagnation:
				print ("\t >>> Resetting the solution to avoid local optima.")
				parent_seq = ''
				parent_eval = []
				randomMutationFlag = True
				stagnation = 0
				maxStagnation *= 1.2

	print ()
	print ()
#	print ("\t Elapsed time:       " + str(time.time() - start) + " secs.")
#	print ("\t ViennaRNA time:     " + str(totalRuntimeVienna) + " secs.")
#	print ("\t Stop:               " + str(stop) + " secs.") 
	return sequences_found, n_evals


############################################################################################################################################################################################################################################
#       MAIN
############################################################################################################################################################################################################################################

start_time = time.time()
n_evals = 0

# Create an instance of the ArgumentParser class
parser = argparse.ArgumentParser(description='Local Search RNA Inverse Folding')

# Add arguments with flags
parser.add_argument('-t', '--target', help='Target structure in dot-bracket notation')
parser.add_argument('-n', '--nSolutions', type=int, default=1, help='Number of RNA sequences (solutions) to be found by the algorithm')
parser.add_argument('--stopping', type=int, default=60,  help='Stopping criterion in seconds')
parser.add_argument('--seed', type=int, help='Seed for random numbers')
parser.add_argument('--stagnationLimit', type=int, default=50, help='Number of iterations with no improvement before resetting the current solution to avoid local optima (by default 50 iterations).')
parser.add_argument('--pm', type=float, default=0.01,  help='Mutation probability (by default 1%)')
parser.add_argument('--maxGC', type=float, default=0.48,  help='Constraint the percentage of GC content (by default 0.48 (48%), as suggested by Kiriakidou et al. [10.1016/j.cell.2007.05.016])')
parser.add_argument('--turner1999', action='store_true', help='Energy Model TURNER1999')
parser.add_argument('--turner2004', action='store_true', help='Energy Model TURNER2004')

# Parse the command line arguments
args = parser.parse_args()

decomposition = False

# Access the values of the arguments
if args.target:
	target = args.target
if args.stopping:
	stop = int(args.stopping)
if args.seed:
	random.seed(int(args.seed))
if args.stagnationLimit:
	stagnationLimit = int(args.stagnationLimit)
if args.nSolutions:
	maxN = int(args.nSolutions)
if args.turner1999:
	# TURNER1999
	p = Path(__file__).resolve().parent
	RNA.params_load(str(p.joinpath("rna_turner1999.par")))
	model = "TURNER1999"
if args.turner2004:
	# TURNER2004
	p = Path(__file__).resolve().parent
	RNA.params_load(str(p.joinpath("rna_turner2004.par")))
	model = "TURNER2004"
if args.maxGC:
	maxGC = float(args.maxGC)
if args.pm:
	pm = float(args.pm)

print ()
print ("-----------------------------------------------------------------------")
print("   _____ _____  ______ ______ _____             _____  _   _           ")
print("  / ____|  __ \|  ____|  ____|  __ \           |  __ \| \ | |   /\     ")
print(" | |  __| |__) | |__  | |__  | |  | |  ______  | |__) |  \| |  /  \    ")
print(" | | |_ |  _  /|  __| |  __| | |  | | |______| |  _  /| . ` | / /\ \   ")
print(" | |__| | | \ \| |____| |____| |__| |          | | \ \| |\  |/ ____ \  ")
print("  \_____|_|  \_\______|______|_____/           |_|  \_\_| \_/_/    \_\ ")
print("       ")
print ("                                       -- A. Rubio-Largo et al. 2023")
print ("-----------------------------------------------------------------------")
print ()
print ("------------------------------------------------------------------------------------------------------------------")
print ("CONFIGURATION: ")
print ()
print (" * TARGET:           "  + target)
print (" * STOPPING:         "  + str(stop) + " secs.")
print (" * SEED:             "  + str(args.seed) )
print (" * NO. OF SOLUTIONS: "  + str(maxN))
print (" * max. GC-cont.:    "  + str(maxGC))
print (" * MODEL:            "  + model)
print (" * MUTATION PROB:    "  + str(pm))
print (" * STAGNATION LIMIT: "  + str(stagnationLimit))
print ("------------------------------------------------------------------------------------------------------------------")

results, n_evals = greed(target, maxGC, maxN, pm, stagnationLimit, stop)

print ()
print ()
print ("RNA sequences found with GC-content <= " + str(maxGC) + ":")
other_rnaS = []
for r in results:
	pGC = calculate_pGC(r,target)
	if (pGC <= maxGC):
		print ("## " + r + " " + str(pGC))
	else:
		other_rnaS.append("@@ " + r + " " + str(pGC))
print()
print ("RNA sequences found with GC-content > " + str(maxGC) + ":")
for r in other_rnaS:
	print (r)

print()
print()
print ("Elapsed time:                   " + str(time.time() - start_time) + " secs.")
print ("Number of function evaluations: " + str(n_evals))
print ('----------------------------------------------------------------------')

print()
print()
print ("---------------------------------------")
print (" FINISHED SUCCESSFULLY ")
print ("---------------------------------------")
sys.exit()
