#!/bin/bash

#
#	Note that the following command line arguments are used specifically for a run of
#	the simulation with relative fitness
#

#0 for relative; 1 for absolute

#General variables
timeSteps=20000
initialPopsize=20000
mud=2.1
chromosomesize=200
numberofchromosomes=23
bentodelratio=0
sb=1
#0 for point; 1 for exponential; and 2 for uniform
bendist=1
#0 for root sb; 1 for single; 2 for root Ncrit
typeofrun=1
#0 for no tskit; 1 for tskit on; 2 for tskit on after burnin
tskitstatus=2
SdtoSbratio=0.029
#0 for Kim et al., 1 for exponential, 2 for point
deldist=1

#
#	The command line arguments below are then specifically used for a run of the simulation
#	with absolute fitness
#

slope=0
seed=24
#rawdata file size in datapoints
rawdatafilesize=10
#change in carrying capacity, type in the difference in popsize
#status of fixation calculation; 0 for OFF; 1 for ON

#
#	=====================================================================================
#	MODULAR EPISTASIS - NOT SUPPORTED IN THIS BUILD. LEAVE BOTH AT 0.
#	=====================================================================================
#	The modular branch of RecombineChromosomesIntoGamete was incomplete (it only ever
#	copied the first half of each chromosome), its elementsperlb indexing overran the
#	gamete buffers, and it has no notion of the modifier-locus mask. The program now
#	REJECTS modularepis=1 at startup rather than silently producing wrong genotypes.
#
#	Both variables must stay defined: they occupy positional argument slots 18 and 19,
#	and removing them would shift every argument after them.
#
#	Original lines, preserved:
#	#0 for runs without modular epistasis; 1 for runs with modular epistasis
#	modularepis=0
#	elementsperl=0
#	=====================================================================================
#

#must remain 0 - modular epistasis is rejected by the program
modularepis=0
#unused when modularepis=0; kept only to fill its positional argument slot
elementsperl=0

#
#	=====================================================================================
#	EVOLUTION OF THE MUTATION RATE - MODIFIER LOCI
#	=====================================================================================
#	Every linkage block carries a fitness component. A fraction of the linkage blocks
#	ALSO carries a "modifier locus", which is in one of two states:
#
#	        +1  mutator
#	        -1  anti-mutator
#
#	A linkage block WITHOUT a modifier locus is a "non-modifier" block: it is permanently
#	held at 0 and can never change state.
#
#	Both mutation rates of an individual are then
#
#	        mu_deleterious = mud * f^n
#	        mu_beneficial  = mub * f^n
#
#	where f = mutator_strength_factor and n is the sum of the modifier states over the
#	whole diploid genome. Modifier states are copied through the SAME recombination
#	breakpoints as the fitness effects, so a modifier allele never separates from its
#	linkage block.
#
#	TO TURN MUTATION-RATE EVOLUTION OFF COMPLETELY: set mutator_switch_rate=0.0 and
#	mutator_strength_factor=1.0 (or simply nmodifierlociperchromosome=0).
#	=====================================================================================
#

#f in mu = mu0 * f^n; 1.0 means mutator alleles have no effect on mutation rate
mutator_strength_factor=1.0
#per-MODIFIER-LOCUS, per-gamete probability of a mutator/anti-mutator switch; 0 disables mutator evolution
#(non-modifier blocks are never tested, so this is now a per-modifier-locus rate, not a per-block rate)
mutator_switch_rate=0.0
#bias applied to the anti-mutator -> mutator switch rate relative to mutator -> anti-mutator
#  anti-mutator -> mutator rate = mutator_switch_rate * mutator_bias
#  mutator -> anti-mutator rate = mutator_switch_rate
#  >1 favours mutators, <1 favours anti-mutators, 1.0 is unbiased
mutator_bias=1.0

#
#	--- Which linkage blocks carry a modifier locus ------------------------------------
#

#Number of modifier loci ON EACH CHROMOSOME. Exactly this many distinct blocks are drawn
#at random within every chromosome, so the total per haploid genome is
#  nmodifierlociperchromosome * numberofchromosomes
#and the modifier loci are spread evenly across the linkage groups. Must be between 0 and
#chromosomesize. 0 means no block carries a modifier locus, so the mutation rate can never
#change. The set is fixed for the whole run and identical in every individual and on both
#homologs; there is no per-individual mask.
nmodifierlociperchromosome=10

#q: fraction of modifier loci that start in the +1 (mutator) state at generation 0.
#  One set of starting positions is drawn and applied to every haplotype, so the founding
#  population is monomorphic at every modifier locus and all later variation is generated
#  by the simulation itself. Every other modifier locus starts at -1 (anti-mutator).
initial_mutator_fraction=0.0

#
#	--- Optional detailed per-individual output ----------------------------------------
#	Writes one row per individual (Wi, log Wi, both realised mutation rates, mutator allele
#	count, anti-mutator count, modifier slots, net modifier sum) into individualtrackingfor*.txt.
#	This is EXPENSIVE: each firing writes $initialPopsize rows. Use trackinterval and
#	trackstartgen to restrict it to the window you actually want to plot.
#	The per-generation POPULATION means and variances are always written to the raw data
#	file regardless of this setting, so leave this off unless you need individual detail.
#

#0 = off, 1 = on. When on, writes individualtracking.txt with one row per individual per
#firing: Generation, Individual, Wi, LogFitness (absolute, offset-independent),
#LogFitnessOffset, both realised mutation rates, MutatorAlleleCount,
#AntiMutatorAlleleCount, ModifierSlots, NetModifierSum, and one MutatorsOnChrN column per
#chromosome. Each firing writes $initialPopsize rows.
trackindividuals=0
#dump every this many N-timesteps (generations). Ignored when trackindividuals=0.
trackinterval=100
#first generation (1-based) eligible for dumping; use this to skip the burn-in
trackstartgen=1
#
#	--- Restart checkpoints -------------------------------------------------------------
#	A checkpoint is a COMPLETE binary image of the simulation: every individual's fitness
#	and modifier-state arrays, the shared modifier-locus positions, the log-fitness offset,
#	the running sum and the burn-in detector's history. It is written to
#	checkpoint_gen<N>.bin inside the output directory.
#
#	SIZE: popsize * 2 * numberofchromosomes * chromosomesize * 12 bytes.
#	  popsize 500,   9200 blocks  ->  about   55 MB per checkpoint
#	  popsize 20000, 9200 blocks  ->  about  2.2 GB per checkpoint
#	Keep this interval MUCH larger than trackinterval, which is far cheaper.
#	0 disables checkpointing entirely.
#
#	TO RESUME: set the MUTATIONLOAD_RESUME environment variable to a checkpoint path, e.g.
#	    MUTATIONLOAD_RESUME=/path/to/checkpoint_gen5000.bin ./mutationload <the 26 arguments>
#	It is an environment variable and not a 27th argument so that this script's argument
#	list never has to change. The checkpoint's popsize, genome length, chromosomesize and
#	numberofchromosomes are checked against the current run and any mismatch aborts.
#
#	CAVEAT: a resumed run is NOT the same realisation as an uninterrupted one. Both random
#	number generators are reseeded deterministically from the stored seed and generation
#	rather than having their internal state restored, so the continuation is reproducible
#	and statistically valid but not bit-identical. Fine for continuing a long run; not
#	suitable for reproducing one exact trajectory.
#
checkpointinterval=0

# This build only runs relative fitness - the fitnesstype argument has been
# removed from the command line - so the string is fixed. Original block:
#if [ $fitnesstype -eq 0 ]
#then
#	fitnessstring="relative_"
#elif [ $fitnesstype -eq 1 ]
#then
#	fitnessstring="absolute_"
#fi
fitnessstring="relative_"

if [ $bendist -eq 0 ]
then
	bendiststring="point_"
elif [ $bendist -eq 1 ]
then
	bendiststring="exponential_"
elif [ $bendist -eq 2 ]
then
	bendiststring="uniform_"
fi

if [ $tskitstatus -eq 0 ]
then
	tskitstatusstring="OFF_"
elif [ $tskitstatus -eq 1 ]
then
	tskitstatusstring="ON_"
elif [ $tskitstatus -eq 2 ]
then
	tskitstatusstring="ON_AFTER_BURNIN_"
fi

if [ $deldist -eq 0 ]
then
	deldiststring="kim_"
elif [ $deldist -eq 1 ]
then
	deldiststring="exponential_"
elif [ $deldist -eq 2 ]
then
	deldiststring="point_"
fi

#mub is written as a formated double in mutation load program
mub=$(echo "$mud * $bentodelratio" | bc -l)
mub=$(printf "%.4f" $mub)

#creates 2 strings; directory refers to the folder where data for the specified parameters will be stored; file is the snapshot of the simulation at its end.
if [ $modularepis -eq 0 ]
then
	directory="datafor_"$fitnessstring"tskitstatus_"$tskitstatusstring"r_"$r"_i_init"$i_init"_s_"$s"_K_"$K"_deldist_"$deldiststring"bendist_"$bendiststring"_mub_"$mub"_chromnum_"$numberofchromosomes"_N0_"$initialPopsize"_mud_"$mud"_L_"$chromosomesize"seed_"$seed"/"
elif [ $modularepis -eq 1 ]
then
	directory="datafor_"$fitnessstring"tskitstatus_"$tskitstatusstring"elementsperlb_"$elementsperl"_r_"$r"_i_init"$i_init"_s_"$s"_K_"$K"_deldist_"$deldiststring"_bendist_"$bendiststring"_mub_"$mub"_chromnum_"$numberofchromosomes"_N0_"$initialPopsize"_mud_"$mud"_L_"$chromosomesize"seed_"$seed"/"
fi

printf "directory path is %s \n" "$directory"

file1="popsnapshotfor_popsize_"$initialPopsize"_tskitstatus_"$tskitstatusstring"mub_"$mub".txt"

#checks if a previous snapshot of the simulation exist. Snapshots are saved as compressed files (.gz) to save space
prevsim=$([ -f $directory$file1".gz" ] && echo 1 || echo 0)

if [ $prevsim -eq 0 ]
then
	snapshot=0
elif [ $prevsim -eq 1 ]
then
	gzip -d $directory$file1".gz"
	snapshot=1
fi


SECONDS=0
echo "start of mutationload program"

# run mutationload program with arguments
./mutationload $timeSteps $initialPopsize $mud $chromosomesize $numberofchromosomes $bentodelratio $sb $bendist $typeofrun $slope $seed $tskitstatus $modularepis $elementsperl $SdtoSbratio $deldist $rawdatafilesize $mutator_strength_factor $mutator_switch_rate $mutator_bias $nmodifierlociperchromosome $initial_mutator_fraction $trackindividuals $trackinterval $trackstartgen $checkpointinterval

echo $SECONDS

echo "end of mutationload program"

#$snapshot $directory$file1

gzip $directory$file1