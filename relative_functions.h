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
double RunSimulationRel(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * Sbname, char * mutator_switch_ratename, char * mutator_biasname, char * mutator_strength_factorname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize, MutatorConfig mutatorconfig, TrackingConfig trackingconfig);

void PerformOneTimeStepRel(int tskitstatus, bool isabsolute, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, int popsize, int totaltimesteps, double currenttimestep, long double *wholepopulationwistree, Individual *wholepopulation, long double * psumofwis, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *parent1gameteFitness, int *parent1gameteMutators, GameteState *parent1state, double *parent2gameteFitness, int *parent2gameteMutators, GameteState *parent2state, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, MutatorConfig mutatorconfig);

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

/* --- output helpers ------------------------------------------------------- */
void WritePopulationModifierSummary(FILE *rawdatafilepointer, Individual *wholepopulation, int popsize, int totalindividualgenomelength, const int *modifierlocuspositions, int nmodifierloci, int *locusmutatorcounts);
void WriteIndividualSnapshot(FILE *individualfilepointer, Individual *wholepopulation, int popsize, int generation, int nmodifierslots);

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
