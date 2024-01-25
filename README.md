# GREED-RNA
### Overview
GREED-RNA is a RNA design algorithm based on a simple greedy evolutionary strategy. It allows to obtain a desired number of sequences solution to the RNA inverse folding problem that meet a given maximum GC-content for the base-pairs, and also returns all the sequences found that solve the structure but do not meet that limit. It also allows to specify a stopping-criterion (maximum execution time), a stagnation limit and a mutation factor.

The main feature is the use of several objective functions (Base-pair distance, Hamming distance, probability over ensemble, partition function, ensemble defect and GC-content of base-pairs) to select the best solution in each iteration, changing their weight according to the problem solving conditions.
In addition, it incorporates a greedy initialization and mutation, which ease obtaining a sequence that folds as the target but also, to expand the search space, the mutation mode changes to random when stagnation is detected.

### Algorithm
Briefly, the process of GREED-RNA follows the next loop:

1. _Greedy initialization_, when the parent sequence is empty.
2. _Mutation_. Default mode is greedy mutation. The position to which it will be applied is randomly selected, but in the case of base pairs only GC or CG are chosen, and an A is always placed at unpaired positions. The mode changes to random when stagnation is detected or the state of the algorithm reaches case 3.
3. _Dynamic selection_. The current sequences of parent and offspring are sorted depending on the state of the algorithm and their values in the objective functions. Both sequences are first sorted according to their value in base-pair distance, thus minimizing it. To break the ties, the other objective functions are used sequentially. The second one is always Hamming distance, but the order of the remaining ones changes depending on the state of the algorithm:
   1. At the beginning of the process: −1·(Probability over ensemble), partition function, ensemble defect, GC-content.
   2. Once a sequence that folds exactly like the target structure has been found (base pair distance = 0): GC-content, partition function, ensemble defect, −1·(probability over ensemble).
   3. When a sequence with a GC-content for the base-pairs lower or equal than the required maximum (maxGC) is found that folds exactly like the target structure, and if the requested number of solutions has not yet been reached: |GC content − maxGC|, partition function, ensemble defect, −1·(probability over ensemble).
      
   The candidate sequence (current parent or offspring) that is located in the first position after sorting will be selected as the parent sequence for the next generation. In this way, the most appropriate criterion is used at all times, instead of always applying the same one.
4. _Stagnation detection_. If after sorting the best (first) candidate sequence is the current parent sequence, the stagnation counter is incremented. When this counter reaches the maximum set, the parent sequence will be reinitialized to avoid local optima, and random mutation will be set. At this point it is also checked if the number of valid solutions found is at least the number of solutions requested. If so, the loop is terminated.

![GREED-RNA-image](https://github.com/iARN-unex/GREED-RNA/assets/118007536/f3fefe23-8836-42b4-8748-003639d5932e)

### Comparative Study.
The performance of GREED-RNA  was evaluated with the widely used Eterna100 benchmark in its V1 and V2 version with their corresponding Turner1999(T99) and Turner2004(T04) energy parameters sets.

For each Eterna100 structure in both versions of Eterna100 GREED-RNA was run 10 times. The parameter settings were: Stopping criterion of 86400 seconds (24 hours), maximum GC-content for the base-pairs of 0.48, mutation probability of 0.01, stagnation limit of 50 and number of solutions of 1. As explained above, the number of solutions applies only to those that meet the maximum GC-content for the base-pairs requirement. Resulting sequences can be found in [data/output](data/output)     [data/V1-T99/Sequences-found](data/V1-T99/Sequences-found) and [data/V2-T04/Sequences-found](data/V2-T04/Sequences-found) folders. Sequences beginning with ## comply with the GC-content for the base-pairs limit, while those beginning with @@ do not.

Results were compared against other RNA inverse folding methods.

### Use guide
To execute GREED-RNA in your local machine the [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/) is needed.

usage: ```GREED-RNA_.py [-h] [-t TARGET] [-n NSOLUTIONS] [--stopping STOPPING] [--seed SEED] [--stagnationLimit STAGNATIONLIMIT] [--pm PM]
                     [--maxGC MAXGC] [--turner1999] [--turner2004]```

Being the parameters:
                     
```
-t, --target. Target structure in dot-bracket notation.
-n, --nSolutions, type=int, default=1. Number of RNA sequences (solutions) to be found by the algorithm.
--stopping, type=int, default=60. Stopping criterion in seconds.
--seed, type=int. Seed for random numbers.
--stagnationLimit, type=int, default=50. Number of iterations with no improvement before resetting the current solution to avoid local optima (by default 50 iterations).
--pm, type=float, default=0.01. Mutation probability (by default 1%).
--maxGC, type=float, default=0.48. Constraint the percentage of GC content (by default 0.48 (48%).
--turner1999. Energy Model TURNER1999.
--turner2004. Energy Model TURNER2004.
```

Only -t/--target and one between --turner1999 or --turner2004 are mandatory.




