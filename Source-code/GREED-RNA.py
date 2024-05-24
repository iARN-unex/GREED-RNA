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
    up_ = []
    U_ = []
    B_ = []
    stack = []
    pair_dict = {}

    for i, char in enumerate(_struct):
        if char == '.':
            up_.append(init + i)
        elif char == '(':
            stack.append(i)
        elif char == ')':
            if not stack:
                raise ValueError('_structure error: Numbers of "(" and ")" do not match')
            opening_index = stack.pop()
            pair_dict[opening_index] = i
        else:
            raise ValueError('_structure error: Wrong character')

    if stack:
        raise ValueError('_structure error: Numbers of "(" and ")" do not match')

    for opening_index in sorted(pair_dict.keys()):
        B_.append(init + pair_dict[opening_index])
        U_.append(init + opening_index)

    U_.extend(up_)

    return B_, U_
    

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
def calculate_pGC_(sequence, structure):
    B_, U_ = getComponents_B_U(structure)

    # Initialize count for CG and GC pairs
    count_CG_GC = 0

    # Loop through the base pairs
    for i in range(len(B_)):
        pair = sequence[U_[i]] + sequence[B_[i]]
        if pair in ('CG', 'GC'):
            count_CG_GC += 1

    # Calculate the pGC value
    total_pairs = len(B_)
    pGC = count_CG_GC / max(total_pairs, 1)

    return pGC

def calculate_pGC(sequence):
    gcCount = sequence.count('G') + sequence.count('C')
    return gcCount / len(sequence)

def calculate_dist_pGC(sequence, minGC, maxGC):
    gcContent = calculate_pGC(sequence)
    return max(minGC - gcContent, 0, gcContent - maxGC)

def hammingDistance(str1, str2):
    if len(str1) != len(str2):
        raise ValueError("Strings must be of the same length")
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

def calculateCompositionM2dRNAs(sequence, structure):
    B_, U_ = getComponents_B_U(structure)

    # Base pairs distribution
    pair_counts = {'AU': 0, 'UA': 0, 'CG': 0, 'GC': 0, 'GU': 0, 'UG': 0}
    for i in range(len(B_)):
        pair = sequence[U_[i]] + sequence[B_[i]]
        if pair in pair_counts:
            pair_counts[pair] += 1

    total_pairs = max(len(B_), 1)
    pAU = (pair_counts['AU'] + pair_counts['UA']) / total_pairs
    pGC = (pair_counts['CG'] + pair_counts['GC']) / total_pairs
    pUG = (pair_counts['GU'] + pair_counts['UG']) / total_pairs
    paired = max(pAU, pGC, pUG)

    # Unpaired bases distribution
    subseq = "".join(sequence[U_[i]] for i in range(len(B_), len(U_)))
    if subseq:
        uA = subseq.count('A') / len(subseq)
        uU = subseq.count('U') / len(subseq)
        uG = subseq.count('G') / len(subseq)
        uC = subseq.count('C') / len(subseq)
        unpaired = max(uA, uU, uG, uC)
    else:
        unpaired = 0

    # Total bases distribution
    total_bases = len(sequence)
    a = sequence.count('A') / total_bases
    u = sequence.count('U') / total_bases
    g = sequence.count('G') / total_bases
    c = sequence.count('C') / total_bases
    total = max(a, u, g, c)

    return round((total + unpaired + paired) * 100, 3), pGC

def evaluateM2DRNAS(s, t):
    # Create a fold compound with the given RNA sequence
    fc = RNA.fold_compound(s)
    
    # Calculate minimum free energy (mfe) and its structure
    mfe_struct, mfe = fc.mfe()
    
    # Calculate base pair distance between the mfe structure and target structure
    f_bpdistance = RNA.bp_distance(mfe_struct, t)
    
    # Rescale the partition function parameters
    fc.exp_params_rescale(mfe)
    
    # Calculate the partition function and its free energy
    pp, f_pf = fc.pf()
    
    # Calculate the ensemble diversity
    f_ensemble = fc.mean_bp_distance()
    
    # Calculate composition and pGC of the mfe structure
    f_composition, pGC = calculateCompositionM2dRNAs(s, mfe_struct)
    
    return f_pf, f_ensemble, f_composition, f_bpdistance, s, mfe_struct, pGC

def evaluate(s, t, minGC, maxGC):
    # Crear un compuesto de plegamiento con la secuencia de ARN dada
    fc = RNA.fold_compound(s)

    # Calcular la energía libre mínima (mfe) y su estructura
    mfe_struct, mfe = fc.mfe()

    # Calcular la distancia de pares de bases entre la estructura mfe y la estructura objetivo
    f_bpdistance = RNA.bp_distance(mfe_struct, t)

    # Reescalar los parámetros de la función de partición
    fc.exp_params_rescale(mfe)

    # Calcular la función de partición y su energía libre
    pp, f_pf = fc.pf()

    # Calcular el defecto del conjunto respecto a la estructura objetivo
    f_ensemble = fc.ensemble_defect(t)

    # Calcular la probabilidad de la estructura objetivo
    f_prob = fc.pr_structure(t)

    # Calcular la distancia de Hamming entre la estructura mfe y la estructura objetivo
    f_dist = hammingDistance(mfe_struct, t)

    # Calcular el valor de pGC para la estructura mfe y estimar su desviación del rango (minGC,maxGC)
    f_pGC = calculate_dist_pGC(s, minGC, maxGC)

    # Devolver los resultados como una tupla, sin redondear
    return mfe_struct, f_bpdistance / len(t), (-1.0) * f_prob, f_pf, f_ensemble, f_pGC, f_dist / len(t)

############################################################################################################################################################################################################################################
#       GREEDY INITIALIZATION
############################################################################################################################################################################################################################################
def greedyInitialization(target, B, U):
    s = ['_'] * len(target)
    triplet = random.choice(['GCA', 'CGA'])
    
    # Fill the base pairs
    for i in range(len(B)):
        if s[B[i]] == '_' and s[U[i]] == '_':
            s[B[i]] = triplet[0]
            s[U[i]] = triplet[1]
    
    # Fill the unpaired positions
    for i in range(len(U)):
        if s[U[i]] == '_':
            s[U[i]] = triplet[2]
    
    return "".join(s)

############################################################################################################################################################################################################################################
#       RANDOM MUTATION
############################################################################################################################################################################################################################################
def randomMutation (sequence, pm, B, U):
	s = list(sequence)
	len_U = len(U)
	len_B = len(B)
	for x in range (pm):
		pos = int(random.uniform(0.00, len_U))
		#pos = random.randint(0, len_U - 1)
		if pos < len_B:
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
	len_U = len(U)
	len_B = len(B)
	for x in range (pm):
		pos = int(random.uniform(0.00, len_U))
		#pos = random.randint(0, len_U - 1)
		if pos < len_B:
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
def greed(target, minGC, maxGC, maxN, pm, maxStagnation, stagnationInc, stop):
    start_time = time.time()
    B, U = getComponents_B_U(target)
    #totalRuntimeVienna = 0

    amount = max(round(len(target) * pm), 1)
    sequences_found = []
    best_sequences_soFar = []
    n_evals = 0
    stagnation = 0
    flagInit = False
    flagBPD = False
    flagGC = False
    randomMutationFlag = False
    minBPD = 999999
    minN = 0

    parent_seq = ''
    parent_eval = []
    iter = 0

    print("\n\n (1) Let's minimize the base-pair distance between predicted and target structures")
    print("    -> Min(Base-pair distance), min(Hamming distance), max(probability over ensemble), min(partition function), min(ensemble defect), min(%GC-content)\n\n")

    while time.time() - start_time < stop:
        iter += 1

        # Greedy initialization (first time and after re-start)
        if parent_seq == '':
            flagInit = True
            if len(best_sequences_soFar) > 0 and random.uniform(0.00, 1.00) > 0.5:
                print("\t\t >>> Resetting the solution from BEST POOL")
                parent_seq = random.choice(best_sequences_soFar)[0]
            else:
                print("\t\t >>> Resetting the solution from GREEDY INITIALIZATION")
                parent_seq = greedyInitialization(target, B, U)
            
            #_initEvalTime = time.time()
            parent_eval = evaluate(parent_seq, target, minGC, maxGC)
            #totalRuntimeVienna += (time.time() - _initEvalTime)
            n_evals += 1

        # Mutation greedy-random
        if not randomMutationFlag:
            offspr_seq = greedyMutation(parent_seq, amount, B, U)
        else:
            offspr_seq = randomMutation(parent_seq, amount, B, U)

        # Evaluate offspring
        if offspr_seq != parent_seq:
            #_initEvalTime = time.time()
            offspr_eval = evaluate(offspr_seq, target, minGC, maxGC)
            #totalRuntimeVienna += (time.time() - _initEvalTime)
            n_evals += 1
        else:
            offspr_eval = parent_eval

        results = [[parent_seq] + list(parent_eval), [offspr_seq] + list(offspr_eval)]
        best_sequences_soFar.append([parent_seq] + list(parent_eval))
        best_sequences_soFar.append([offspr_seq] + list(offspr_eval))

        # Remove duplicates and sort
	# SELECTION DYNAMICALLY
	#  (0) seq
	#  (1) mfe_struct
	#  (2) f_bpdistance/len(t)
	#  (3) f_prob
	#  (4) f_p
	#  (5) f_ensemble
	#  (6) f_pGC (distance to range minGC, maxGC)
	#  (7) f_dist/len(t)
        if flagBPD:
            if not flagGC:
                results.sort(key=lambda x: (x[2], x[7], x[6], x[4], x[5], x[3]))
                best_sequences_soFar = sorted(
                    [list(x) for x in set(tuple(x) for x in best_sequences_soFar)],
                    key=lambda x: (x[2], x[7], x[6], x[4], x[5], x[3])
                )[:50]
            else:
                results.sort(key=lambda x: (x[2], x[7], abs(x[6]-maxGC), x[4], x[5], x[3]))
                best_sequences_soFar = sorted(
                    [list(x) for x in set(tuple(x) for x in best_sequences_soFar)],
                    key=lambda x: (x[2], x[7], abs(x[6]-maxGC), x[4], x[5], x[3])
                )[:50]
        else:
            results.sort(key=lambda x: (x[2], x[7], x[3], x[5], x[4], x[6]))
            best_sequences_soFar = sorted(
                [list(x) for x in set(tuple(x) for x in best_sequences_soFar)],
                key=lambda x: (x[2], x[7], x[3], x[5], x[4], x[6])
            )[:50]

        if parent_seq != results[0][0] or flagInit:
            if parent_eval[1] > results[0][2] or (parent_eval[1] == 0 and parent_eval[5] > results[0][6]) or flagInit:
                print('\t ->' + str([parent_seq] + list(parent_eval)))

            flagInit = False
            stagnation = 0
            parent_seq = results[0][0]
            parent_eval = results[0][1:]

            if parent_eval[1] == 0:
                if not flagBPD:
                    flagBPD = True
                    print("\n\n    -> Found sequence! Elapsed time: " + str(time.time() - start_time) + " secs.\n\n")
                    print(" (2) We found a sequence that folds exactly like the target structure. Let's minimize the GC-content")
                    print("    -> ORDER: Min(Base-pair distance), min(Hamming distance), min(GC-content), min(partition function), min(ensemble defect), max(probability over ensemble)\n\n")

                if parent_seq not in sequences_found:
                    sequences_found.append(parent_seq)

                if parent_eval[5] == 0.0: # pGC is in the range (minGC, maxGC):
                    flagGC = True
                    randomMutationFlag = True
                    minN += 1
                    if minN < maxN:
                        print("\n\n    -> Satisfied GC-constraint! Elapsed time (2): " + str(time.time() - start_time) + " secs.\n\n")
                        print(" (3) We found 1 sequence that folds exactly like the target structure with a GC-content in the range [" + str(minGC) + ", " + str(maxGC) + "]. Let's find " + str(maxN - minN) + " more structures")
                        print("    -> ORDER: Min(Base-pair distance), min(Hamming distance), min(abs(GC-content-maxGC)), min(partition function), min(ensemble defect), max(probability over ensemble)\n\n")

            if minBPD > parent_eval[1]:
                minBPD = parent_eval[1]
        else:
            stagnation += 1
            if flagBPD and flagGC and minN >= maxN:
                print("\n\n    -> Finished!! Elapsed time: " + str(time.time() - start_time) + " secs.")
                break
            elif stagnation >= maxStagnation:
                print("\t >>> Resetting the solution to avoid local optima. (" + str(maxStagnation) + " * " + str(stagnationInc) + " = " + str(maxStagnation * stagnationInc) + ")")
                parent_seq = ''
                parent_eval = []
                randomMutationFlag = True
                stagnation = 0
                maxStagnation *= stagnationInc

    return sequences_found, n_evals

############################################################################################################################################################################################################################################
#       MAIN
############################################################################################################################################################################################################################################

def main():
    start_time = time.time()
    n_evals = 0

    # Create an instance of the ArgumentParser class
    parser = argparse.ArgumentParser(description='Local Search RNA Inverse Folding')

    # Add arguments with flags
    parser.add_argument('-t', '--target', required=True, help='Target structure in dot-bracket notation')
    parser.add_argument('-n', '--nSolutions', type=int, default=1, help='Number of RNA sequences (solutions) to be found by the algorithm')
    parser.add_argument('--stopping', type=int, default=60, help='Stopping criterion in seconds')
    parser.add_argument('--seed', type=int, help='Seed for random numbers')
    parser.add_argument('--stagnationLimit', type=int, default=50, help='Number of iterations with no improvement before resetting the current solution to avoid local optima (by default 50 iterations).')
    parser.add_argument('--pm', type=float, default=0.01, help='Mutation probability (by default 1%)')
    parser.add_argument('--maxGC', type=float, default=0.6, help='Constraint the percentage of GC content (by default min. 60%)')
    parser.add_argument('--minGC', type=float, default=0.4, help='Constraint the percentage of GC content (by default min. 40%)')
    parser.add_argument('--stagnationInc', type=float, default=0.20, help='Constraint the increase of iterations for the next number of iterations with no improvements (stagnation)')
    parser.add_argument('--turner1999', action='store_true', help='Energy Model TURNER1999')
    parser.add_argument('--turner2004', action='store_true', help='Energy Model TURNER2004')

    # Parse the command line arguments
    args = parser.parse_args()

    # Access the values of the arguments
    target = args.target
    stop = int(args.stopping)
    if args.seed:
        random.seed(int(args.seed))
    stagnationLimit = int(args.stagnationLimit)
    maxN = int(args.nSolutions)
    maxGC = float(args.maxGC)
    minGC = float(args.minGC)
    if maxGC < minGC:
      print ("*** WARNING *** maxGC < minGC; setting maxGC to minGC")
      maxGC = minGC

    pm = float(args.pm)
    stagnationInc = 1 + float(args.stagnationInc)

    if args.turner1999:
        # TURNER1999
        p = Path(__file__).resolve().parent
        RNA.params_load(str(p.joinpath("rna_turner1999.par")))
        model = "TURNER1999"
    elif args.turner2004:
        # TURNER2004
        p = Path(__file__).resolve().parent
        RNA.params_load(str(p.joinpath("rna_turner2004.par")))
        model = "TURNER2004"
    else:
        model = "Default"

    print()
    print("-----------------------------------------------------------------------")
    print("   _____ _____  ______ ______ _____             _____  _   _           ")
    print("  / ____|  __ \|  ____|  ____|  __ \           |  __ \| \ | |   /\     ")
    print(" | |  __| |__) | |__  | |__  | |  | |  ______  | |__) |  \| |  /  \    ")
    print(" | | |_ |  _  /|  __| |  __| | |  | | |______| |  _  /| . ` | / /\ \   ")
    print(" | |__| | | \ \| |____| |____| |__| |          | | \ \| |\  |/ ____ \  ")
    print("  \_____|_|  \_\______|______|_____/           |_|  \_\_| \_/_/    \_\ ")
    print("       ")
    print("                                       -- A. Rubio-Largo et al. 2023")
    print("-----------------------------------------------------------------------")
    print()
    print("------------------------------------------------------------------------------------------------------------------")
    print("CONFIGURATION: ")
    print()
    print(f" * TARGET:           {target}")
    print(f" * STOPPING:         {stop} secs.")
    print(f" * SEED:             {args.seed}")
    print(f" * NO. OF SOLUTIONS: {maxN}")
    print(f" * max. GC-cont.:    {maxGC}")
    print(f" * min. GC-cont.:    {minGC}")
    print(f" * MODEL:            {model}")
    print(f" * MUTATION PROB:    {pm}")
    print(f" * STAGNATION LIMIT: {stagnationLimit}")
    print(f" * STAGNATION INC.:  {args.stagnationInc}")
    print("------------------------------------------------------------------------------------------------------------------")

    results, n_evals = greed(target, minGC, maxGC, maxN, pm, stagnationLimit, stagnationInc, stop)

    print()
    print()
    print(f"RNA sequences found with GC-content in the range [{minGC},{maxGC}]:")
    other_rnaS = []
    for r in results:
        pGC = calculate_pGC(r)
        if pGC >= minGC and pGC <= maxGC:
            print(f"## {r} {pGC}")
        else:
            other_rnaS.append(f"@@ {r} {pGC}")
    print()
    print(f"RNA sequences found with GC-content lower than {minGC} or higher than {maxGC}:")
    for r in other_rnaS:
        print(r)

    print()
    print()
    print(f"Elapsed time:                   {time.time() - start_time} secs.")
    print(f"Number of function evaluations: {n_evals}")
    print('----------------------------------------------------------------------')
    print()
    print()
    print("---------------------------------------")
    print(" FINISHED SUCCESSFULLY ")
    print("---------------------------------------")
    sys.exit()

if __name__ == "__main__":
    main()
