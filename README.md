# HMMCoriolan
HMM Coding of Coriolanus for final project.

First part is coding of 2-block states of Coriolanus, a transition matrix and calculating theoretical/experimental coding of a concatenated state arithmetic code.

Second part is a Hidden Markov Model with the same states, and emission probabilities, using hmmtrain to refine the transmission/emission matrices. This, however, ran into a bunch of issues, including:
1. Long convergence of hmmtrain
2. Refined transmission matrix implying that all but one of the states loop back to itself infinitely, and the leading eigenvalue not being a valid probability vector.

Update: Well... there are less issues. Although the experimental encoding is much worse then the theoretical.
