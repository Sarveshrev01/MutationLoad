#ifndef RELATIVE_FUNCTIONS_H_INCLUDED
#define RELATIVE_FUNCTIONS_H_INCLUDED 1

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
#include "sharedfunc_flag.h"
#include "main.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

/* MutatorConfig (see sharedfunc_flag.h) now carries the modifier-locus settings:
 * f, the switch rate, the bias, the NUMBER of modifier loci per chromosome and
 * the initial mutator fraction. TrackingConfig carries the per-individual dump
 * settings. Both are passed by value. */
double RunSimulationRel(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * Sbname, char * mutator_switch_ratename, char * mutator_biasname, char * mutator_strength_factorname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize, int randomnumberseed, MutatorConfig mutatorconfig, TrackingConfig trackingconfig);

void PerformOneTimeStepRel(int tskitstatus, bool isabsolute, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, int popsize, int totaltimesteps, double currenttimestep, long double *wholepopulationwistree, Individual *wholepopulation, long double * psumofwis, long double logfitnessoffset, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *parent1gameteFitness, int *parent1gameteMutators, GameteState *parent1state, double *parent2gameteFitness, int *parent2gameteMutators, GameteState *parent2state, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, MutatorConfig mutatorconfig);

void InitializePopulationRel(int tskitstatus, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, long double *wholepopulationwistree, Individual *wholepopulation, int popsize, int totalpopulationgenomelength, int chromosomesize, int numberofchromosomes, int totaltimesteps, long double * psumofwis, int *modifierlocuspositions, int *pnmodifierloci, MutatorConfig mutatorconfig, FILE *miscfilepointer);

int ChooseVictim(int populationsize);
int ChooseParentWithTree(long double *wholepopulationwistree, int popsize, long double sumofwis, FILE *miscfilepointer);

/* --- modifier-locus helpers (mutation-rate evolution) ---------------------
 * DrawModifierLocusPositions draws EXACTLY mutatorconfig.lociperchromosome
 * distinct blocks on EACH chromosome, so the modifier loci are spread evenly
 * across the linkage groups. The positions are haploid indices, returned in
 * ascending order, and are shared by every individual and both homologs for the
 * whole run. Returns the total number drawn.
 * ------------------------------------------------------------------------- */
int DrawModifierLocusPositions(int *positions, int chromosomesize, int numberofchromosomes, int lociperchromosome);
void SeedInitialMutatorStates(int *stateshaploid, const int *positions, int nmodifierloci, int haploidgenomelength, double initialmutatorfraction);

/* --- tree-sequence bootstrap (item 7: honour tskitstatus == 2) ------------ */
void SeedTreeSequenceTables(tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, int popsize, int haploidgenomelength, double nodetime);

/* -------------------------------------------------------------------------
 * RENORMALISATION (fitness rescaling)
 * -------------------------------------------------------------------------
 * Wi = exp(logFitness - logfitnessoffset). Because parents are drawn with
 * probability proportional to Wi, multiplying every Wi by a common constant
 * leaves every selection probability exactly unchanged. Subtracting a constant
 * from every individual's log-fitness does precisely that, so the offset can be
 * moved at any time without perturbing the dynamics.
 *
 * RenormalizeFitness re-centres the offset on the population's current mean
 * log-fitness, recomputes every Wi, and rebuilds the Fenwick selection tree and
 * sumofwis from scratch. Called at most once per generation, and only when the
 * mean has drifted past RENORMALIZATION_THRESHOLD.
 * ------------------------------------------------------------------------- */
void RenormalizeFitness(Individual *wholepopulation, int popsize, long double newoffset, long double *wholepopulationwistree, long double *psumofwis, double mutator_strength_factor, double baseline_deleterious_rate, double baseline_beneficial_rate);

/* Population-wide log-fitness summary, computed in a single O(N) pass:
 * the mean, the maximum and the minimum of the ABSOLUTE logFitness. */
void SummariseLogFitness(Individual *wholepopulation, int popsize, long double *pmean, long double *pmax, long double *pmin);

/* -------------------------------------------------------------------------
 * RESTART CHECKPOINTS
 * -------------------------------------------------------------------------
 * A checkpoint is a complete image of the simulation state: every individual's
 * fitness and modifier-state arrays, the shared modifier-locus positions, the
 * log-fitness offset, the running sum, the burn-in detector's history, and all
 * the parameters needed to validate a resume.
 *
 * SIZE. The body is popsize * 2L * 12 bytes. At popsize 500 with a 9200-block
 * genome that is about 55 MB per checkpoint; at popsize 20000 it is about
 * 2.2 GB. Set checkpointinterval accordingly - it is deliberately separate from
 * the per-individual tracking interval, which is far cheaper.
 *
 * RESUMING. Set the environment variable MUTATIONLOAD_RESUME to a checkpoint
 * path. This is an environment variable rather than a command-line argument so
 * that the 26-argument layout, and therefore every submission script, stays
 * unchanged. The header is validated against the current parameters and the run
 * aborts on any mismatch.
 *
 * A RESUMED RUN IS NOT THE SAME REALISATION as an uninterrupted one. Both
 * generators are reseeded deterministically from the stored seed and generation
 * rather than having their internal state restored, so the continuation is
 * reproducible and statistically valid but not bit-identical. That is fine for
 * continuing a long run; it is not suitable for reproducing one exact trajectory.
 *
 * Checkpoints are binary and are only guaranteed readable by the same build on
 * the same architecture.
 * ------------------------------------------------------------------------- */
int WriteCheckpoint(const char *path, Individual *wholepopulation, int popsize, int totalindividualgenomelength, int chromosomesize, int numberofchromosomes, int Nxtimesteps, int generation, int randomnumberseed, const int *modifierlocuspositions, int nmodifierloci, int isburninphaseover, int endofburninphase, int endofdelay, int Nxtimestepsafterburnin, double currenttimestep, double deleteriousmutationrate, double beneficialmutationrate, double Sb, double Sd, MutatorConfig mutatorconfig, long double logfitnessoffset, long double sumofwis, const double *literallyjustlast200Ntimesteps, const double *last200Ntimestepsvariance, const double *logaveragefitnesseachNtimesteps, FILE *miscfilepointer);

int ReadCheckpoint(const char *path, Individual *wholepopulation, int popsize, int totalindividualgenomelength, int chromosomesize, int numberofchromosomes, int *pgeneration, int *prandomnumberseed, int *modifierlocuspositions, int *pnmodifierloci, int *pisburninphaseover, int *pendofburninphase, int *pendofdelay, int *pNxtimestepsafterburnin, double *pcurrenttimestep, long double *plogfitnessoffset, long double *psumofwis, double *literallyjustlast200Ntimesteps, double *last200Ntimestepsvariance, double *logaveragefitnesseachNtimesteps, FILE *miscfilepointer);

/* --- output helpers ------------------------------------------------------- */
void WritePopulationModifierSummary(FILE *rawdatafilepointer, Individual *wholepopulation, int popsize, int totalindividualgenomelength, const int *modifierlocuspositions, int nmodifierloci, int *locusmutatorcounts);
void WriteIndividualSnapshot(FILE *individualfilepointer, Individual *wholepopulation, int popsize, int generation, int nmodifierslots, long double logfitnessoffset, const int *modifierlocuspositions, int nmodifierloci, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, int *chromosomemutatorcounts);

/* ---------------------------------------------------------------------------
 * CalculateWi - COMMENTED OUT (item 10: completely unused).
 * ---------------------------------------------------------------------------
 * Superseded by the inline Wi recomputation inside PerformBirth() and by
 * UpdateIndividual(), both of which work on the Individual struct rather than on
 * a pair of raw gamete arrays. Kept here, commented, rather than deleted.
 *
 * double CalculateWi(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength);
 * ------------------------------------------------------------------------- */

#endif // RELATIVE_FUNCTIONS_H_INCLUDED
