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
#include "relative_functions.h"
/* #include "absolute_functions.h"  - COMMENTED OUT (Tier 3): absolute_functions.c
 * is not in the Makefile and nothing in this file references any symbol it
 * declares. Including it only pulled in declarations for a translation unit that
 * is never linked. */
#include "sharedfunc_flag.h"
#include "main.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

/* -------------------------------------------------------------------------
 * MAXMUTATIONSPERGAMETE  (fix for the fixed-size Sds[] buffer)
 * -------------------------------------------------------------------------
 * The deleterious effect sizes drawn for one gamete are held in a fixed-size
 * stack buffer so that no allocation happens on the hot path. The old buffer
 * was 30 entries, which was safe only because the mutation rate was constant.
 * Now that mu evolves upwards (mu = mu0 * f^n) an unlucky individual can draw
 * far more than 30 mutations, so the buffer is raised to 10000 entries
 * (10000 * sizeof(double) = 80 kB of stack, well inside the default 8 MB).
 *
 * If a draw ever exceeds this, the run ABORTS with a diagnostic rather than
 * silently truncating the mutations, because a truncated draw would bias the
 * realised mutation rate without any visible symptom. Hitting this limit means
 * mutator_strength_factor is too large for the chosen baseline rate.
 * ------------------------------------------------------------------------- */
#define MAXMUTATIONSPERGAMETE 10000

void MutateGamete(int tskitstatus, int isburninphaseover,  tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationsitesarray, tsk_id_t childnode, int totaltimesteps, double currenttimestep, bool isabsolute, int totalindividualgenomelength, double *gamete, GameteState *gs, double mutationeffectsize)
{
    /* NOTE (Tier 3): the isburninphaseover parameter is unused in this function.
     * In the original code it gated recording for absolute runs only; the relative
     * path now controls recording through the flag passed in the tskitstatus slot
     * (see item 7). The parameter is kept so this shared signature is unchanged. */
    tsk_id_t idofnewmutation;
    int mutatedsite = DetermineMutationSite(totalindividualgenomelength/2);
    double appliedeffect;
    if(isabsolute){
        appliedeffect = mutationeffectsize;
    }else{
        appliedeffect = log(1 + mutationeffectsize);
    }
    gamete[mutatedsite] += appliedeffect;
    /* Keep the gamete's running log-fitness in step, so PerformBirth never has
     * to re-sum the genome. */
    if (gs != NULL) gs->logfitnesssum += (long double) appliedeffect;
    char derivedstate[400];
    sprintf(derivedstate, "%.11f", mutationeffectsize);
    /* NOTE (item 7): callers in the relative path now pass a *recording-active*
     * flag in the tskitstatus slot rather than the raw tskitstatus, so that
     * tskitstatus == 2 ("record only after burn-in") is honoured here without
     * changing this function's contract. See PerformOneTimeStepRel(). */
    if (tskitstatus != 0){
        idofnewmutation = tsk_mutation_table_add_row(&treesequencetablecollection->mutations, wholepopulationsitesarray[mutatedsite], childnode, TSK_NULL, ((double) totaltimesteps - currenttimestep), derivedstate, 12, NULL, 0);
        check_tsk_error(idofnewmutation); 
    }
}

// PerformDeath remains mostly unchanged, just fixing array references passed as NULL in relative runs
double PerformDeath(bool isabsolute, int tskitstatus, int isburninphaseover, int maxPopSize, int *pPopSize, int victim, int deleteriousdistribution, long double *wholepopulationselectiontree, Individual *wholepopulation, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double b_0, double r,  int i_init, double s, long double *psumofload, long double *psumofloadsquared, tsk_id_t * wholepopulationnodesarray, FILE *miscfilepointer)
{
    /* The isabsolute branch was an empty stub. main() aborts on absolute runs and
     * RunSimulationRel refuses them, so isabsolute is always false here and only
     * the relative body ever executed. Original structure preserved:
     *
     * if(isabsolute){
     *     // Absolute logic
     * }
     * else{
     */
    *psumofloads -= wholepopulation[victim].fitness;
    wholepopulation[victim].fitness = 0.0;
    /* } */
    Fen_set(wholepopulationselectiontree, maxPopSize, 0.0, victim);
    return 0.0;
}

/* NOTE (Tier 1): ismodular and elementsperlb are unused - modular epistasis is
 * not supported in this build and main() rejects it. The parameters are kept so
 * this shared signature is unchanged. */
void PerformBirth(int tskitstatus, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t childnode1, tsk_id_t childnode2, bool isabsolute, double *parent1gameteFitness, int *parent1gameteMutators, const GameteState *parent1state, double *parent2gameteFitness, int *parent2gameteMutators, const GameteState *parent2state, int maxPopSize, int *pPopSize, int birthplace, Individual *wholepopulation, int totalindividualgenomelength, int deleteriousdistribution, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double b_0, double r,  int i_init, double s, long double *psumofload, long double *psumofloadsquared, FILE *miscfilepointer, double mutator_strength_factor, double baseline_deleterious_rate, double baseline_beneficial_rate)
{
    int i;
    long double newwi;
    int halfgenome = totalindividualgenomelength/2;

    // Copy gametes into the individual at 'birthplace'
    for (i = 0; i < halfgenome; i++) {
        wholepopulation[birthplace].fitnessArray[i] = parent1gameteFitness[i];
        wholepopulation[birthplace].fitnessArray[halfgenome + i] = parent2gameteFitness[i];
        wholepopulation[birthplace].mutatorArray[i] = parent1gameteMutators[i];
        wholepopulation[birthplace].mutatorArray[halfgenome + i] = parent2gameteMutators[i];
    }

    /* -------------------------------------------------------------------
     * O(1) assembly of the newborn's summary state.
     * -------------------------------------------------------------------
     * Both gametes already carry their own running log-fitness and modifier
     * sum, built while RecombineChromosomesIntoGamete copied them and updated
     * by MutateGamete and SwitchModifierLoci. The diploid values are simply
     * the two halves added together, so neither this function nor the rate
     * update needs to sweep the genome. The two O(2L) sweeps that used to live
     * here and in UpdateIndividual are gone.
     * ------------------------------------------------------------------- */
    wholepopulation[birthplace].logFitness     = parent1state->logfitnesssum + parent2state->logfitnesssum;
    wholepopulation[birthplace].netModifierSum = parent1state->modifiersum   + parent2state->modifiersum;

    /* The isabsolute branch was an empty stub. main() aborts on absolute runs and
     * RunSimulationRel refuses them, so isabsolute is always false here and only
     * the relative body ever executed. Original structure preserved:
     *
     * if(isabsolute){
     *     // Absolute logic...
     * }
     * else{
     */
    newwi = expl(wholepopulation[birthplace].logFitness);

    Fen_set(wholepopulationselectiontree, maxPopSize, newwi, birthplace);
    wholepopulation[birthplace].fitness = (double) newwi;
    *psumofloads += newwi;
    /* } */

    // Update Cached Mutation Rate for new individual. O(1): reads only the
    // already-assembled logFitness and netModifierSum.
    RefreshIndividualRates(&wholepopulation[birthplace], mutator_strength_factor, baseline_deleterious_rate, baseline_beneficial_rate);

    /* NOTE (item 7): as in MutateGamete, the relative path passes a
     * recording-active flag here, not the raw tskitstatus. */
    if (tskitstatus != 0){
        wholepopulationnodesarray[birthplace*2] = childnode1;
        wholepopulationnodesarray[birthplace*2 + 1] = childnode2; 
    }
}

Individual createIndividual(double *fitnessArray, int *mutatorArray, int totalindividualgenomelength){
    Individual ind;
    ind.fitnessArray = fitnessArray;
    ind.mutatorArray = mutatorArray;
    ind.logFitness = 0.0;
    ind.fitness = 1.0;
    ind.mutationRate = 0.0;
    ind.beneficialMutationRate = 0.0;
    ind.netModifierSum = 0;
    return ind;
}

/* -------------------------------------------------------------------------
 * RefreshIndividualRates   -   O(1)
 * -------------------------------------------------------------------------
 * Derives everything that depends on the two maintained summary values:
 *
 *      fitness        = exp(logFitness)
 *      mu_deleterious = mu_d0 * f^n          <- item 8: both rates use the same
 *      mu_beneficial  = mu_b0 * f^n             f and the same n
 *
 * where n = netModifierSum, the NET sum over the diploid genome (+1 per mutator,
 * -1 per anti-mutator, 0 per non-modifier block).
 *
 * This replaces the old UpdateIndividual, which swept all 2L blocks on every
 * single birth. Both logFitness and netModifierSum are now carried through
 * recombination and mutation incrementally, so nothing here touches the arrays.
 * ------------------------------------------------------------------------- */
void RefreshIndividualRates(Individual *ind, double mutator_strength_factor, double baseline_deleterious_rate, double baseline_beneficial_rate){
    double modifierfactor = pow(mutator_strength_factor, (double) ind->netModifierSum);
    ind->fitness = (double) expl(ind->logFitness);
    ind->mutationRate           = baseline_deleterious_rate * modifierfactor;
    ind->beneficialMutationRate = baseline_beneficial_rate  * modifierfactor;
}

/* -------------------------------------------------------------------------
 * RecomputeIndividualFromArrays   -   O(2L)
 * -------------------------------------------------------------------------
 * Rebuilds logFitness and netModifierSum from the arrays themselves, then
 * refreshes the derived values. Used when initialising the founding population,
 * and available as the exact re-sync for the incrementally maintained values if
 * floating-point drift over a very long run ever needs correcting.
 * ------------------------------------------------------------------------- */
void RecomputeIndividualFromArrays(Individual *ind, int totalindividualgenomelength, double mutator_strength_factor, double baseline_deleterious_rate, double baseline_beneficial_rate){
    int i;
    long double sum = 0.0;
    int net = 0;
    for (i = 0; i < totalindividualgenomelength; i++){
        sum += ind->fitnessArray[i];
        net += ind->mutatorArray[i];   /* +1, -1 or 0 - no mask lookup needed */
    }
    ind->logFitness = sum;
    ind->netModifierSum = net;
    RefreshIndividualRates(ind, mutator_strength_factor, baseline_deleterious_rate, baseline_beneficial_rate);
}

void RecombineChromosomesIntoGamete(bool isabsolute, int tskitstatus, bool ismodular, int elementsperlb, int isburninphaseover, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int persontorecombine, int chromosomesize, int numberofchromosomes, double *gameteFitness, int *gameteMutators, GameteState *gs, Individual *wholepopulation, int totalindividualgenomelength)
{
    int recombinationsite, startchromosome, h, i, returnvaluefortskit;
    
    tsk_id_t parentnode1 = (tsk_id_t) 2*persontorecombine;
    tsk_id_t parentnode2 = (tsk_id_t) (2*persontorecombine + 1);
    int chromatid_len = totalindividualgenomelength / 2;

    /* The gamete's summary state is rebuilt from scratch for every gamete. */
    gs->nmutatorpositions = 0;
    gs->nantimutatorpositions = 0;
    gs->logfitnesssum = 0.0;
    gs->modifiersum = 0;

    /* NOTE (item 7): the relative path passes a recording-active flag here. */
    if (tskitstatus != 0){
        *childnode = tsk_node_table_add_row(&treesequencetablecollection->nodes, 0, ((double) totaltimesteps - currenttimestep), TSK_NULL, TSK_NULL, NULL, 0);
        check_tsk_error(*childnode);
    }

    for (h = 0; h < numberofchromosomes; h++) {
        startchromosome = pcg32_boundedrand(2); 
        do {
            recombinationsite = pcg32_boundedrand(chromosomesize);
        } while (recombinationsite == 0); 
        
        if (tskitstatus != 0){
            if (startchromosome == 0){
                returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (double)(h*chromosomesize), (double)(h*chromosomesize + recombinationsite), wholepopulationnodesarray[parentnode1], *childnode, NULL, 0);
                check_tsk_error(returnvaluefortskit);
                returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (h*chromosomesize + recombinationsite), ((h+1)*chromosomesize), wholepopulationnodesarray[parentnode2], *childnode, NULL, 0);
                check_tsk_error(returnvaluefortskit);
            }else{
                returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (h*chromosomesize), (h*chromosomesize + recombinationsite), wholepopulationnodesarray[parentnode2], *childnode, NULL, 0);
                check_tsk_error(returnvaluefortskit);
                returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (h*chromosomesize + recombinationsite), ((h+1)*chromosomesize), wholepopulationnodesarray[parentnode1], *childnode, NULL, 0);
                check_tsk_error(returnvaluefortskit);
            }
        }
        
        // Recombine Fitness AND Mutator Arrays
        if(!ismodular){
            /* The fitness effect, the modifier state and (in inherited-mask mode)
             * the modifier mask of a linkage block are all copied through the SAME
             * breakpoint, so a modifier allele can never be separated from its
             * linkage block. While we are already walking every position we also
             * build the modifier index used by the switching step - this is what
             * makes the switching step O(number of modifier loci) instead of
             * O(genome length). See ProduceMutatedGamete()/SwitchModifierLoci(). */
            for (i = 0; i < chromosomesize; i++) {
                int source_offset;
                int idx = h*chromosomesize + i;
                int state;

                if (i < recombinationsite) {
                    source_offset = (startchromosome == 0) ? 0 : chromatid_len;
                } else {
                    source_offset = (startchromosome == 0) ? chromatid_len : 0;
                }

                gameteFitness[idx]  = wholepopulation[persontorecombine].fitnessArray[source_offset + idx];
                state               = wholepopulation[persontorecombine].mutatorArray[source_offset + idx];
                gameteMutators[idx] = state;

                /* Running log-fitness, so PerformBirth never re-sums the genome. */
                gs->logfitnesssum += (long double) gameteFitness[idx];

                /* The state array is self-describing: non-zero means this block
                 * carries a modifier locus. +1 mutator, -1 anti-mutator, 0 inert.
                 * No mask array is consulted or carried. */
                if (state == 1) {
                    gs->mutatorpositions[gs->nmutatorpositions++] = idx;
                    gs->modifiersum += 1;
                } else if (state == -1) {
                    gs->antimutatorpositions[gs->nantimutatorpositions++] = idx;
                    gs->modifiersum -= 1;
                }
            }
        } else {
             /* ---------------------------------------------------------------
              * MODULAR-EPISTASIS BRANCH - DISABLED.
              * ---------------------------------------------------------------
              * This project does not use modular epistasis (run with
              * modularepis = 0). The original block below was incomplete: it
              * only ever copied the FIRST half of each chromosome ("... similar
              * for second half" was never written) and its elementsperlb
              * indexing overran the gamete buffers. It is also not aware of the
              * modifier mask. Rather than leave code that would silently produce
              * wrong genotypes, the branch now aborts. The original lines are
              * preserved verbatim, commented out, immediately below.
              * ---------------------------------------------------------------
             for (i = 0; i < recombinationsite*elementsperlb; i++) {
                int source_offset = (startchromosome == 0) ? 0 : chromatid_len;
                int idx = h*chromosomesize*elementsperlb + i;
                gameteFitness[idx] = wholepopulation[persontorecombine].fitnessArray[source_offset + idx];
                // Assuming mutators align with elements per lb, or are just per block?
                // Based on context, mutators seem to be per linkage block.
                // If elementsperlb > 1, the mutator array indexing needs careful handling or mutators need to be per element.
                // Assuming 1-to-1 mapping for simplicity given the provided code context.
                gameteMutators[idx] = wholepopulation[persontorecombine].mutatorArray[source_offset + idx];
            }
            // ... similar for second half
              * --------------------------------------------------------------- */
            fprintf(stderr, "Error: modular epistasis (modularepis = 1) is not supported by the mutation-rate-evolution build. Run with modularepis = 0.\n");
            exit(1);
        }
    }
}

/* -------------------------------------------------------------------------
 * SwitchModifierLoci   (fix for the per-locus Bernoulli loop)
 * -------------------------------------------------------------------------
 * WHAT CHANGED AND WHY
 *
 * The old implementation looped over EVERY position of the gamete (~4600) and
 * drew one uniform random number per position to decide whether that position
 * switched state. For a 20000 x 20000 run that is roughly 3.7e12 random draws
 * spent almost entirely on positions that cannot switch at all.
 *
 * Two changes, neither of which alters the model:
 *
 *   1. Only modifier loci are considered. Non-modifier blocks are held at 0 and
 *      can never change state, so testing them was pure waste.
 *      RecombineChromosomesIntoGamete already walks the gamete to copy it, so it
 *      builds the two position lists (+1 and -1) for free.
 *
 *   2. Instead of one Bernoulli trial per eligible locus, the NUMBER of
 *      switches is drawn once from the exact Binomial distribution that those
 *      independent trials define, and then that many DISTINCT loci are chosen
 *      uniformly at random via a partial Fisher-Yates shuffle. This is
 *      distributionally identical to the per-locus loop but costs
 *      O(number of switches) instead of O(number of loci).
 *
 * Both directions use the rates the original code used:
 *      anti-mutator -> mutator :  mutator_switch_rate * mutator_bias
 *      mutator -> anti-mutator :  mutator_switch_rate
 *
 * Both draws are taken against the PRE-switch state, exactly as the original
 * per-locus loop did (each locus was evaluated once, against the state it had
 * on entry), so a locus can never be flipped twice in one call.
 *
 * Anti-mutators are always stored as -1, so each flip moves the gamete's
 * modifier sum by exactly +2 or -2 and gs->modifiersum is kept in step.
 * ------------------------------------------------------------------------- */
static void SwitchModifierLoci(int *gameteMutators, GameteState *gs, double mutator_switch_rate, double mutator_bias, gsl_rng * randomnumbergeneratorforgamma)
{
    unsigned int numberofswitches;
    unsigned int j;
    double uprate, downrate;

    if (gs == NULL || mutator_switch_rate <= 0.0) return;

    uprate   = mutator_switch_rate * mutator_bias;
    downrate = mutator_switch_rate;
    if (uprate   > 1.0) uprate   = 1.0;
    if (downrate > 1.0) downrate = 1.0;

    /* anti-mutator (-1) -> mutator (+1). Each flip moves n by +2. */
    if (gs->nantimutatorpositions > 0 && uprate > 0.0) {
        numberofswitches = gsl_ran_binomial(randomnumbergeneratorforgamma, uprate, (unsigned int) gs->nantimutatorpositions);
        for (j = 0; j < numberofswitches; j++) {
            /* Partial Fisher-Yates: swap a uniformly chosen not-yet-used entry
             * into slot j, guaranteeing j distinct positions after j steps. */
            int remaining = gs->nantimutatorpositions - (int) j;
            int pick = (int) j + (int) pcg32_boundedrand((uint32_t) remaining);
            int chosen = gs->antimutatorpositions[pick];
            gs->antimutatorpositions[pick] = gs->antimutatorpositions[j];
            gs->antimutatorpositions[j] = chosen;
            gameteMutators[chosen] = 1;
            gs->modifiersum += 2;
        }
    }

    /* mutator (+1) -> anti-mutator (-1). Each flip moves n by -2. */
    if (gs->nmutatorpositions > 0 && downrate > 0.0) {
        numberofswitches = gsl_ran_binomial(randomnumbergeneratorforgamma, downrate, (unsigned int) gs->nmutatorpositions);
        for (j = 0; j < numberofswitches; j++) {
            int remaining = gs->nmutatorpositions - (int) j;
            int pick = (int) j + (int) pcg32_boundedrand((uint32_t) remaining);
            int chosen = gs->mutatorpositions[pick];
            gs->mutatorpositions[pick] = gs->mutatorpositions[j];
            gs->mutatorpositions[j] = chosen;
            gameteMutators[chosen] = -1;
            gs->modifiersum -= 2;
        }
    }
}

bool ProduceMutatedGamete(int tskitstatus, int isburninphaseover, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int parent, bool isabsolute, int individualgenomelength, double parent_specific_deleterious_rate, double parent_specific_beneficial_rate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *gameteFitness, int *gameteMutators, GameteState *gs, double mutator_switch_rate, double mutator_bias, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer)
{
    int k, numberofbeneficialmutations, numberofdeleteriousmutations;
    double generatedSb;
    /* See MAXMUTATIONSPERGAMETE at the top of this file: raised from 30 to 10000
     * because the mutation rate now evolves. */
    static double Sds[MAXMUTATIONSPERGAMETE];

    // 1. Fitness Mutations (Uses Parent's Mutation Rate)
    bool stayInWhileLoop = true;
    while (stayInWhileLoop) {
        stayInWhileLoop = false;
        numberofdeleteriousmutations = DetermineNumberOfMutations(parent_specific_deleterious_rate);

        /* Hard stop rather than a silent truncation - a truncated draw would
         * bias the realised mutation rate with no visible symptom. */
        if (numberofdeleteriousmutations > MAXMUTATIONSPERGAMETE) {
            fprintf(miscfilepointer, "\nFATAL: %d deleterious mutations drawn for one gamete of individual %d, which exceeds MAXMUTATIONSPERGAMETE (%d).\n", numberofdeleteriousmutations, parent, MAXMUTATIONSPERGAMETE);
            fprintf(miscfilepointer, "The parent's realised deleterious mutation rate was %g. Reduce mutator_strength_factor or the baseline mutation rate, or raise MAXMUTATIONSPERGAMETE in sharedfunc_flag.c.\n", parent_specific_deleterious_rate);
            fflush(miscfilepointer);
            fprintf(stderr, "FATAL: deleterious mutation count %d exceeds MAXMUTATIONSPERGAMETE (%d). See miscellaneous.txt.\n", numberofdeleteriousmutations, MAXMUTATIONSPERGAMETE);
            exit(1);
        }

        for (k = 0; k < numberofdeleteriousmutations; k++) {
            if (deleteriousdistribution == 0) {
                 Sds[k] = (gsl_ran_gamma(randomnumbergeneratorforgamma, 0.169, 1327.4)/23646);
            } else if (deleteriousdistribution == 1) {
                 Sds[k] = gsl_ran_exponential(randomnumbergeneratorforgamma, Sd);
            } else {
                 Sds[k] = Sd;
            }
            if (!isabsolute && Sds[k] >= 1) {
                stayInWhileLoop = true;
                break;
            }
        }
    }

    for (k = 0; k < numberofdeleteriousmutations; k++) {
        double effect = isabsolute ? Sds[k] : -Sds[k];
        MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gameteFitness, gs, effect);
    }
    
    /* ---------------------------------------------------------------------
     * BENEFICIAL MUTATIONS (item 8)
     * ---------------------------------------------------------------------
     * Two fixes here:
     *  (a) the number of beneficial mutations is now drawn from the PARENT'S
     *      realised beneficial rate mu_b0 * f^n, i.e. the same modifier
     *      equation as the deleterious rate with a different baseline, rather
     *      than from the unmodified population-wide beneficialmutationrate; and
     *  (b) beneficialdistribution is honoured again. The previous version had
     *      collapsed every distribution to a point effect of Sb ("Simplified
     *      for brevity"). The three branches below are restored from the
     *      pre-refactor implementation:
     *          0 -> point,       effect = Sb
     *          1 -> exponential, mean  = Sb
     *          2 -> uniform on [0, 2*Sb]
     * Effect sign convention is unchanged: negative under absolute fitness,
     * positive under relative fitness.
     * --------------------------------------------------------------------- */
    numberofbeneficialmutations = DetermineNumberOfMutations(parent_specific_beneficial_rate);

    if (beneficialdistribution == 0) {
        //point distribution
        for (k = 0; k < numberofbeneficialmutations; k++) {
            generatedSb = Sb;
            MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gameteFitness, gs, (isabsolute ? -generatedSb : generatedSb));
        }
    } else if (beneficialdistribution == 1) {
        //exponential distribution
        for (k = 0; k < numberofbeneficialmutations; k++) {
            generatedSb = gsl_ran_exponential(randomnumbergeneratorforgamma, Sb);
            MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gameteFitness, gs, (isabsolute ? -generatedSb : generatedSb));
        }
    } else if (beneficialdistribution == 2) {
        //uniform distribution
        for (k = 0; k < numberofbeneficialmutations; k++) {
            double upperlimitforuniform = (2 * Sb);
            generatedSb = gsl_ran_flat(randomnumbergeneratorforgamma, 0, upperlimitforuniform);
            MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gameteFitness, gs, (isabsolute ? -generatedSb : generatedSb));
        }
    } else {
        fprintf(miscfilepointer, "Error: type of distribution for beneficial effect sizes not recognized.");
        fflush(miscfilepointer);
        exit(0);
    }

    /* 2. Modifier-locus switching (-1 <-> +1).
     *    Only modifier loci are eligible; non-modifier blocks stay at 0 forever.
     *    See the comment block above SwitchModifierLoci for what changed. */
    SwitchModifierLoci(gameteMutators, gs, mutator_switch_rate, mutator_bias, randomnumbergeneratorforgamma);

    return true;
}
