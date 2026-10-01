# Introduction
This project implements a Gibbs sampler to identify short, shared DNA motifs across a set of sequences — a common problem in molecular biology where a regulatory signal (e.g., a transcription factor binding site or ribosome binding site) is known to exist somewhere in each sequence, but its exact position and composition aren't known in advance.

Because real genomic datasets are far too large to exhaustively check every possible motif position across every sequence, we use Gibbs Sampling - a Markov Chain Monte Carlo (MCMC) approach.

The algorithm starts from a random guess at the motif's location in each sequence, then iteratively refines those guesses: on each iteration, one sequence is set aside, a position weight matrix (PWM) is built from the current guesses in every other sequence, and the set-aside sequence's guess is updated by scoring all possible windows against that PWM and sampling a new position probabilistically (rather than always taking the best-scoring window). Repeating this thousands of times allows the guesses to converge on the sequences' shared motif, without ever exhaustively searching the full solution space.

We test our implementation on two datasets: (1) promoter regions upstream of Bacillus subtilis coding sequences, pre-filtered for a fragment of the Shine-Dalgarno motif, and (2) NRF1 ChIP-seq peak sequences, to recover the motif associated with NRF1 transcription factor binding.

# Pseudocode


```


```

# Successes

# Struggles

# Personal Reflections
## Group Leader

## Other members
Dianah: Project 3 is my first time being a collaborator instead of a project leader, so I am learning that side of GitHub as I go, trying to figure out forking, 
how to open a pull request into someone else’s branch instead of my own and what my responsibilities look like when I am not managing the whole repo. 
I am still getting my footing with it, and I think it’s helping me understand GitHub better.

# Generative AI Appendix
