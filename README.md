# GREED-RNA
### Overview
GREED-RNA is a RNA design algorithm based on a simple greedy evolutionary strategy.

The main feature is the use of several objective functions (Base-pair distance, Hamming distance, probability over ensemble, partition function, ensemble defect and GC-content of base-pairs) to select the best solution in each iteration, changing their weight according to the problem solving conditions.
In addition, it incorporates a greedy initialization and mutation, which ease obtaining a sequence that folds as the target but also, to expand the search space, the mutation mode changes to random when stagnation is detected.

### Algorithm
Briefly, the process of GREED-RNA follows the next loop:

1. _Greedy initialization_, when the parent sequence is empty.
2. _Mutation_. Default mode is greedy mutation. The position to which it will be applied is randomly selected, but in the case of base pairs only GC or CG are chosen, and an A is always placed at unpaired positions. The mode changes to random when stagnation is detected or the state of the algorithm reaches case 3.
3. _Dynamic selection_. As explained above, the current sequences of parent and offspring are sorted depending on the state of the algorithm and their values in the objective functions.

![GREED-RNA-image](https://github.com/iARN-unex/GREED-RNA/assets/118007536/f3fefe23-8836-42b4-8748-003639d5932e)
