#ifndef SHAREFUNC_FLAG_H_INCLUDED
#define SHAREFUNC_FLAG_H_INCLUDED 1

#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>
#include <err.h>
#include "dependencies/pcg_basic.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

/* =========================================================================
 * MODIFIER-LOCUS MODEL (mutation-rate evolution)
 * =========================================================================
 * Every linkage block carries a fitness component. A USER-CHOSEN NUMBER of the
 * linkage blocks ON EACH CHROMOSOME additionally carries a "modifier locus".
 * A modifier locus is in one of two states; a block without a modifier locus is
 * a "non-modifier" block and is permanently inert.
 *
 *      state  +1   -> mutator allele
 *      state  -1   -> anti-mutator allele
 *      state   0   -> non-modifier block, NEVER changes
 *
 * The state array is therefore SELF-DESCRIBING: mutatorArray[i] != 0 means
 * position i carries a modifier locus. No separate mask array is needed.
 *
 * The realised mutation rates of an individual are
 *
 *      mu_deleterious = mu_d0 * f^n
 *      mu_beneficial  = mu_b0 * f^n
 *
 * with f = mutator_strength_factor and n = Individual.netModifierSum, the sum of
 * the modifier states over the WHOLE DIPLOID genome. Because anti-mutators are
 * -1, n is a NET SUM: an anti-mutator-heavy individual gets mu below mu_d0.
 * Non-modifier blocks contribute 0 by construction.
 *
 * WHICH blocks carry a modifier locus is fixed for the whole run and identical
 * in every individual and on both homologs. The positions are drawn once at
 * startup (see DrawModifierLocusPositions) and held in a single run-level array.
 *
 * Because the modifier state lives in the same per-block array as the fitness
 * effect and is copied with the SAME recombination breakpoints, a modifier
 * allele can never separate from the linkage block it sits on.
 * ========================================================================= */

/* -------------------------------------------------------------------------
 * MutatorConfig
 * -------------------------------------------------------------------------
 * Everything the user can set about the modifier-locus / mutation-rate-evolution
 * model. Filled once in main() from the command line and then passed by value
 * down to RunSimulationRel().
 * ------------------------------------------------------------------------- */
typedef struct{
    double strengthfactor;         /* f in mu = mu0 * f^n. 1.0 = mutator alleles have no effect. */
    double switchrate;             /* Per-MODIFIER-LOCUS, per-gamete probability of switching     */
                                   /* state. 0 disables modifier evolution entirely.             */
    double bias;                   /* Multiplier applied to the anti-mutator -> mutator rate      */
                                   /* relative to the mutator -> anti-mutator rate.              */
    int    lociperchromosome;      /* Number of modifier loci ON EACH CHROMOSOME. Exactly this    */
                                   /* many distinct blocks are drawn per chromosome, so the total */
                                   /* per haploid genome is lociperchromosome * numberofchromosomes */
                                   /* and they are spread evenly across the linkage groups.       */
    double initialmutatorfraction; /* q: fraction of modifier loci that start in the +1 (mutator)  */
                                   /* state at generation 0. One set of starting positions is      */
                                   /* drawn and applied to every haplotype, so the founding        */
                                   /* population is monomorphic at every modifier locus.          */
} MutatorConfig;

/* -------------------------------------------------------------------------
 * TrackingConfig
 * -------------------------------------------------------------------------
 * Controls the OPTIONAL per-individual dump and (from a later commit) the
 * restart checkpoint.
 * ------------------------------------------------------------------------- */
typedef struct{
    int enabled;            /* 0 = off, 1 = on                                               */
    int interval;           /* Dump every this many N-timesteps (generations). Must be >= 1.  */
    int startgen;           /* First generation (1-based) eligible for dumping.               */
    int checkpointinterval; /* RESERVED, not yet acted on. Will write a restart checkpoint    */
                            /* every this many generations; 0 = never.                       */
} TrackingConfig;

/* Global Individual struct definition */
typedef struct{
    double *fitnessArray; // Still need to change the name to WiArray
    int *mutatorArray;    // Per-block modifier STATE: +1 mutator, -1 anti-mutator,
                          // 0 non-modifier. Non-zero therefore means "is a modifier locus".
    long double logFitness; // Running sum of fitnessArray over the whole diploid genome.
                          // Maintained INCREMENTALLY (see GameteState) so that a birth costs
                          // O(1) instead of a full O(2L) sweep.
    long double fitness;  // Overall fitness (Wi) = expl(logFitness - logfitnessoffset).
                          // long double to match logFitness and sumofwis: Wi is the
                          // exponential of a sum of thousands of block effects, and a
                          // double truncation here was both a precision mismatch with
                          // the Fenwick tree (long double) and an earlier overflow
                          // point than the tree itself. On x86-64 this is the 80-bit
                          // type; on Apple Silicon long double IS double, so the run
                          // behaves identically there.
    double mutationRate;  // Actual DELETERIOUS mutation rate based on modifier loci (mu_d0 * f^n)
    double beneficialMutationRate; // Actual BENEFICIAL mutation rate (mu_b0 * f^n)
    int netModifierSum;   // n, the sum of mutatorArray over the whole diploid genome.
                          // The mutator and anti-mutator COUNTS are not stored: with a fixed,
                          // shared set of modifier loci the number of modifier slots is a
                          // run-level constant M, so
                          //     mutatorCount     = (M + n) / 2
                          //     antimutatorCount = (M - n) / 2
} Individual;

/* -------------------------------------------------------------------------
 * GameteState
 * -------------------------------------------------------------------------
 * Everything one gamete carries besides its two arrays. Built for free by
 * RecombineChromosomesIntoGamete, which already walks every gamete position to
 * copy it, then updated incrementally by MutateGamete and SwitchModifierLoci,
 * and finally consumed by PerformBirth.
 *
 * This is what makes a birth O(1) rather than O(genome length): the old code
 * summed all 2L fitness blocks in PerformBirth and then summed all 2L of them
 * AGAIN in UpdateIndividual. Both sweeps are gone.
 *
 * The two position lists restrict modifier switching to the modifier loci, so
 * the switching step costs O(number of switches) instead of O(genome length).
 * ------------------------------------------------------------------------- */
typedef struct{
    int *mutatorpositions;      /* gamete positions whose state is +1                      */
    int  nmutatorpositions;
    int *antimutatorpositions;  /* gamete positions whose state is -1                      */
    int  nantimutatorpositions;
    long double logfitnesssum;  /* running sum of this gamete's fitness effects            */
    int  modifiersum;           /* running sum of this gamete's modifier states            */
} GameteState;

/* MutateGamete now also accumulates the applied effect into gs->logfitnesssum,
 * so the gamete's running log-fitness stays correct without a re-sweep. */
void MutateGamete(int tskitstatus, int isburninphaseover, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationsitesarray, tsk_id_t childnode, int totaltimesteps, double currenttimestep, bool isabsolute, int totalindividualgenomelength, double *gamete, GameteState *gs, double mutationeffectsize);

double PerformDeath(bool isabsolute, int tskitstatus, int isburninphaseover, int maxPopSize, int *pPopSize, int victim, int deleteriousdistribution, long double *wholepopulationselectiontree, Individual *wholepopulation, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double b_0, double r,  int i_init, double s, long double *psumofload, long double *psumofloadsquared, tsk_id_t * wholepopulationnodesarray, FILE *miscfilepointer);

/* PerformBirth takes the two gametes' GameteState instead of a modifier mask, and
 * assembles the newborn's logFitness and netModifierSum by ADDING them - no sweep. */
void PerformBirth(int tskitstatus, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t childnode1, tsk_id_t childnode2, bool isabsolute, double *parent1gameteFitness, int *parent1gameteMutators, const GameteState *parent1state, double *parent2gameteFitness, int *parent2gameteMutators, const GameteState *parent2state, int maxPopSize, int *pPopSize, int birthplace, Individual *wholepopulation, int totalindividualgenomelength, int deleteriousdistribution, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double b_0, double r,  int i_init, double s, long double *psumofload, long double *psumofloadsquared, FILE *miscfilepointer, long double logfitnessoffset, double mutator_strength_factor, double baseline_deleterious_rate, double baseline_beneficial_rate);

// Individual helper functions
Individual createIndividual(double *fitnessArray, int *mutatorArray, int totalindividualgenomelength);

/* Derives fitness and both mutation rates from the already-maintained logFitness
 * and netModifierSum. O(1); does NOT touch the arrays. */
void RefreshIndividualRates(Individual *ind, long double logfitnessoffset, double mutator_strength_factor, double baseline_deleterious_rate, double baseline_beneficial_rate);

/* Full O(2L) recomputation of logFitness and netModifierSum straight from the
 * arrays. Used at initialisation, and available to re-sync the incrementally
 * maintained values if floating-point drift ever needs correcting. */
void RecomputeIndividualFromArrays(Individual *ind, int totalindividualgenomelength, long double logfitnessoffset, double mutator_strength_factor, double baseline_deleterious_rate, double baseline_beneficial_rate);

void RecombineChromosomesIntoGamete(bool isabsolute, int tskitstatus, bool ismodular, int elementsperlb, int isburninphaseover, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int persontorecombine, int chromosomesize, int numberofchromosomes, double *gameteFitness, int *gameteMutators, GameteState *gs, Individual *wholepopulation, int totalindividualgenomelength);

bool ProduceMutatedGamete(int tskitstatus, int isburninphaseover, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int parent, bool isabsolute, int individualgenomelength, double parent_specific_deleterious_rate, double parent_specific_beneficial_rate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *gameteFitness, int *gameteMutators, GameteState *gs, double mutator_switch_rate, double mutator_bias, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer);

#endif // SHAREFUNC_FLAG_H_INCLUDED
