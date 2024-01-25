# GREED-RNA
### Overview
GREED-RNA is a RNA design algorithm based on a simple greedy evolutionary strategy.

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
4. _Stagnation detection_. If after sorting the best (first) candidate sequence is the current parent sequence, the stagnation counter is incremented. When this counter reaches the maximum set, the parent sequence will be reinitialized to avoid local optima, and random mutation will be set. At this point it is also checked if the number
of valid solutions found is at least the number of solutions requested. If so, the loop is terminated.

![GREED-RNA-image](https://github.com/iARN-unex/GREED-RNA/assets/118007536/f3fefe23-8836-42b4-8748-003639d5932e)
