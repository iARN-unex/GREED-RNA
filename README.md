# GREED-RNA

N. Lozano-García, Á. Rubio-Largo and J. M. Granado-Criado, "A Simple yet Effective Greedy Evolutionary Strategy for RNA Design," in *IEEE Transactions on Evolutionary Computation*, doi: [10.1109/TEVC.2024.3461509](https://doi.org/10.1109/TEVC.2024.3461509).

### Overview
GREED-RNA is a RNA design algorithm based on a simple greedy evolutionary strategy. It allows to obtain a desired number of sequences solution to the RNA inverse folding problem  within a desired range of GC-content, and also returns all the sequences found that solve the structure but do not meet that requirement. It also allows to specify a stopping-criterion (maximum execution time), a stagnation limit and a mutation factor.

The main feature is the use of several objective functions (Base-pair distance, Hamming distance, probability over ensemble, partition function, ensemble defect and GC-content distance) to select the best solution in each iteration, changing their weight according to the problem solving conditions.
In addition, it incorporates a greedy initialization and mutation, which ease obtaining a sequence that folds as the target but also, to expand the search space, the mutation mode changes to random when stagnation is detected.

### Algorithm
Briefly, the process of GREED-RNA follows the next loop:

1. _Greedy initialization_, when the parent sequence is empty. All base-pairs are initialized as GC or CG and unpaired positions are all initialized as A. 
2. _Mutation_. Default mode is greedy mutation. The position to which it will be applied is randomly selected, but in the case of base pairs only GC or CG are chosen, and an A is always placed at unpaired positions. The mode changes to random when stagnation is detected or the state of the algorithm reaches case 3.
3. _Dynamic selection_. The current sequences of parent and offspring are sorted depending on the state of the algorithm and their values in the objective functions. Both sequences are first sorted according to their value in base-pair distance, thus minimizing it. To break the ties, the other objective functions are used sequentially. The second one is always Hamming distance, but the order of the remaining ones changes depending on the state of the algorithm:
   1. At the beginning of the process: −1·(Probability over ensemble), partition function, ensemble defect, GC-content distance.
   2. Once a sequence that folds exactly like the target structure has been found (base pair distance = 0): GC-content distance, partition function, ensemble defect, −1·(probability over ensemble).
   3. When a sequence with a GC-content within the required range is found that folds exactly like the target structure, and if the requested number of solutions has not yet been reached: GC-content distance, partition function, ensemble defect, −1·(probability over ensemble).
      
   The candidate sequence (current parent or offspring) that is located in the first position after sorting will be selected as the parent sequence for the next generation. In this way, the most appropriate criterion is used at all times, instead of always applying the same one.
4. _Stagnation detection_. If after sorting the best (first) candidate sequence is the current parent sequence, the stagnation counter is incremented. When this counter reaches the maximum set, the parent sequence will undergo a restart to prevent being trapped in local optima, and random mutation will be set. The restart randomly chooses between a new greedily initialized solution and one of the Top50 solutions found, with an equal probability for each option. At this point it is also checked if the number of valid solutions found is at least the number of solutions requested. If so, the loop is terminated.

<!-- ![GREED-RNA-image](https://github.com/iARN-unex/GREED-RNA/assets/118007536/f3fefe23-8836-42b4-8748-003639d5932e)-->
![GREED-RNA-image](https://github.com/iARN-unex/GREED-RNA/blob/main/algorithm.png)

### Comparative Study.
The performance of GREED-RNA  was evaluated with the widely used Eterna100 benchmark in its V1 and V2 version with their corresponding Turner1999(T99) and Turner2004(T04) energy parameters sets.

For each Eterna100 structure in both versions of Eterna100 GREED-RNA was run 10 times. The parameter settings were: Stopping criterion of 86400 seconds (24 hours), GC-content range: [40, 60]%, mutation probability of 0.5%, stagnation limit of 50, stagnation increment of 5% and number of solutions of 1. As explained above, the number of solutions applies only to those that meet the GC-content range requirement. Resulting sequences can be found in [data/output](data/output)     [data/V1-T99/Sequences-found](data/V1-T99/Sequences-found) and [data/V2-T04/Sequences-found](data/V2-T04/Sequences-found) folders. Sequences beginning with ## meet the GC-content range requirement, while those beginning with @@ do not.

Results were compared against other RNA inverse folding methods.

### Use guide
To execute GREED-RNA in your local machine the [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/) is needed.

usage: ```GREED-RNA.py [-h] -t TARGET [-n NSOLUTIONS] [--stopping STOPPING] [--seed SEED] [--stagnationLimit STAGNATIONLIMIT] [--pm PM] [--maxGC MAXGC] [--minGC MINGC] [--stagnationInc STAGNATIONINC] [--turner1999] [--turner2004]```
                     


Being the parameters:
                     
```
-t, --target. Target structure in dot-bracket notation.
-n, --nSolutions, type=int, default=1. Number of RNA sequences (solutions) to be found by the algorithm.
--stopping, type=int, default=60. Stopping criterion in seconds.
--seed, type=int. Seed for random numbers.
--stagnationLimit, type=int, default=50. Number of iterations with no improvement before resetting the current solution to avoid local optima (by default 50 iterations).
--pm, type=float, default=0.01. Mutation probability (by default 1%).
--maxGC, type=float, default=0.6. Constraint the percentage of GC content (by default max. 60%).
--minGC, type=float, default=0.4. Constraint the percentage of GC content (by default min. 40%).
--stagnationInc', type=float, default=0.20. Constraint the increase of iterations for the next number of iterations with no improvements (stagnation))
--turner1999. Energy Model TURNER1999.
--turner2004. Energy Model TURNER2004.
```

Only -t/--target and one between --turner1999 or --turner2004 are mandatory.




