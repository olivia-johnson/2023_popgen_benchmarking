import sys
import os
import msprime
import pyslim
import numpy as np
import tskit
import pandas as pd
import time
import itertools
import allel
import statistics
import tracemalloc


## SET PARAMETERS ##
tmpdir="/2023_popgen_benchmarking/"
s_pop=10000  ## scaled population size
sequenceSize=1e6  ## simulated sequence length (1Mb)
recRate=1e-6    ## scaled recombination rate
mutRate=1e-7    ## scaled mutation rate
replicates=10   ## number of replicates

## COALESCENT BURN-IN ##
burnin_time=[]
memory=[]
diversity=[]
simtype=[]
for i in range(replicates):
    # starting the memory monitoring
    tracemalloc.start()
    
    # start time of msprime coalescent simulation
    start_time = time.time()
    
    ## simulate burn-in tree seqeunce
    burnin_ts = msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize,)

    # end time of msprime coalescent simulation
    end_time=time.time()
    # add peak memory consumption to vector IN BYTES
    memory.append(tracemalloc.get_traced_memory()[1])
    tracemalloc.stop()
    
    # add time of tree only simulation and label to vector
    burnin_time.append(end_time-start_time)
    simtype.append("msprime_ts")
    
    # convert tree to SLiM input and output to folder to be read into forward simulation
    burnin_ts = pyslim.annotate(burnin_ts, model_type="WF", tick=1,    stage="late")
    burnin_ts.dump("{0}/burnin_{1}.trees".format(tmpdir,i))
    
    # start memroy monitoring again
    tracemalloc.start()
    # start time for overlaying mutation
    muts_time = time.time()
    
    # overlay neutral mutations on tree
    burnin_ts=msprime.sim_mutations(burnin_ts, rate=mutRate)
    # end of time to overlay mutations
    mute_time=time.time()
    # add peak memory consumption to vector IN BYTES
    memory.append(tracemalloc.get_traced_memory()[1])
    tracemalloc.stop()
    # sum time for generating tree and overlaying mutation to vector, add label to vector
    burnin_time.append((end_time-start_time)+(mute_time-muts_time))
    simtype.append("msprime_mut")
    # calculate diversity and add to vectors
    div=burnin_ts.diversity(mode="site")
    diversity.append(div)
    diversity.append(div)

## combine vectord into table
data = {'time': burnin_time, 'memory': memory, 'diversity': diversity, 'sim_type':simtype}
bench_data = pd.DataFrame(data)
# export table as txt file
bench_data.to_csv("{0}/burnin_benchmarking_msprime.txt".format(tmpdir), index=False, header = True, sep = "\t")

## SLIM CLASSICAL BURN-IN
for i in range(replicates):
    # run SLiM classical simulation
        # time and memory consumption output by SLiM as well as diversity
    os.system("slim -d sim_run={0} ~/2023_popgen_benchmarking/burnin_mut.slim ".format(i))

## SLIM TREE SEQEUNCE RECORDING (10N GENERATIONS) BURN-IN
sts_diversity=[]
for i in range(replicates):
    # run SLiM tree sequence simulation 
    os.system("slim -d sim_run={0} ~/2023_popgen_benchmarking/burnin_ts.slim".format(i))
    
  ## CALCULATE DIVERSITY
    # load tree output from SLiM, format, and overlay neutral mutations mutations
    burnin=tskit.load("{0}burnin_ts_{1}.trees".format(tmpdir,i)).simplify()
    burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1,    stage="late")
    burnin_ts=msprime.sim_mutations(burnin_ts, rate=mutRate)
    # calculate diversity
    div=burnin_ts.diversity(mode="site")
    # append diversity value to vector
    sts_diversity.append(div)
    # read in results from SLiM and add diversity column
sts=pd.read_csv('{0}burnin_benchmarking_slim_ts.txt'.format(tmpdir), sep="\t", header=0)
sts.insert(3, "diversity", sts_diversity)
# export complete table as txt file
sts.to_csv("{0}burnin_benchmarking_slim_ts.txt".format(tmpdir), index=False, header = True, sep = "\t")

## SLIM TREE SEQUENCE USING (checkCoalescence=T) BURN-IN
stscc_diversity=[]
for i in range(replicates):
    # run SLiM simulation withe checkCoalescence = T
    os.system("slim -d sim_run={0} ~/2023_popgen_benchmarking/burnin_tscc.slim".format(i))
    
  ## CALCULATE DIVERSITY
    # load tree output from SLiM, format, and overlay neutral mutations mutations
    burnin=tskit.load("{0}burnin_tscc_{1}.trees".format(tmpdir,i)).simplify()
    burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1,    stage="late")
    burnin_ts=msprime.sim_mutations(burnin_ts, rate=mutRate)
    # calculate diversity
    div=burnin_ts.diversity(mode="site")
    stscc_diversity.append(div)
    # read in results from SLiM and add diversity column
sts=pd.read_csv('{0}burnin_benchmarking_slim_tscc.txt'.format(tmpdir), sep="\t", header=0)
sts.insert(3, "diversity", stscc_diversity)
# export complete table as txt file
sts.to_csv("{0}burnin_benchmarking_slim_tscc.txt".format(tmpdir), index=False, header = True, sep = "\t")



## FORWARD SIMULATIONS WITH SELECTION USING TREE SEQUENCE RECORDING ##

fts_time=[]
fts_sim_type=[]
fts_model=[]
for i in range(replicates):
    # start time of simulation
    start_time = time.time()
    # run SLiM tree sequence simulation with single locus model
    os.system("slim -d sim_run={0} ~/2023_popgen_benchmarking/forward_single_locus.slim".format(str(i)) )
    # load in resulting tree sequence
    slim_ts=tskit.load("{0}/treeseq_single_{1}.trees".format(tmpdir,i)).simplify()
    # overlay neutral mutations
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)
    # end time of simulation
    end_time=time.time()
    # add runtime, simulation type and selection model to vector
    fts_time.append(end_time-start_time)
    fts_sim_type.append("ts")
    fts_model.append("single_locus")
    
for i in range(replicates):
    # start time of simulation
    start_time = time.time()
    # run SLiM tree sequence simulation with multilocus model
    os.system("slim -d sim_run={0} ~/2023_popgen_benchmarking/forward_multi_locus.slim".format(str(i)))
    # load in resulting tree sequence
    slim_ts=tskit.load("{0}/treeseq_multi_{1}.trees".format(tmpdir,i)).simplify()
    # overlay neutral mutations
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)
    # end time of simulation
    end_time=time.time()
    # add runtime, simulation type and selection model to vector
    fts_time.append(end_time-start_time)
    fts_sim_type.append("ts")
    fts_model.append("multilocus")
# create table of resource usage for forward simulations using tree sequence recording
smdata = {'sim_type': fts_sim_type,  'model': fts_model,'time': fts_time}
smbench_data = pd.DataFrame(smdata)
# output table to txt file
smbench_data.to_csv("{0}/forward_benchmarking_time_ts.txt".format(tmpdir), index=False, header = True, sep = "\t")

## CLASSICAL FORWARD SIMULATIONS OF SELECTION
fm_time=[]
fm_sim_type=[]
fm_model=[]
for i in range(replicates):
    # start time of simulation
    start_time = time.time()
    # run classical SLiM simulation with multilocus model
    os.system("slim -d sim_run={0} ~/2023_popgen_benchmarking/forward_multilocus_mut.slim".format(str(i)) )
    # end time of simulation
    end_time=time.time()
    # add runtime, simulation type and selection model to vector
    fm_time.append(end_time-start_time)
    fm_sim_type.append("mut")
    fm_model.append("multilocus")
# create table of resource usage for classical forward simulations of multilocus selection model
smdata = {'sim_type': fm_sim_type,  'model': fm_model,'time': fm_time}
smbench_data = pd.DataFrame(smdata)
# output table to txt file
smbench_data.to_csv("{0}/forward_mut_benchmarking_time_multi.txt".format(tmpdir), index=False, header = True, sep = "\t")

fm_time=[]
fm_sim_type=[]
fm_model=[]
for i in range(replicates):
    # start time of simulation
    start_time = time.time()
    # run classical SLiM simulation with single locus model
    os.system("slim -d sim_run={0} ~/2023_popgen_benchmarking/forward_single_locus_mut.slim".format(str(i)) )
    # end time of simulation
    end_time=time.time()
    # add runtime, simulation type and selection model to vector
    fm_time.append(end_time-start_time)
    fm_sim_type.append("mut")
    fm_model.append("single_locus")
# create table of resource usage for classical forward simulations of single locus selection model
smdata = {'sim_type': fm_sim_type,  'model': fm_model,'time': fm_time}
smbench_data = pd.DataFrame(smdata)
# output table to txt file
smbench_data.to_csv("{0}/forward_mut_benchmarking_time_single.txt".format(tmpdir), index=False, header = True, sep = "\t")


## ANALYSIS OF SIMUALTION DATA ##

calc_time=[]
total_time=[]
fts_stat_type=[]
fts_model=[]
fts_memory=[]
fts_stat=[]
dt=[]
for i in range(replicates):
    # load slingle locus tree sequence data
    slim_ts=tskit.load("{0}/treeseq_single_{1}.trees".format(tmpdir,i)).simplify()
    # overlay neutral mutations
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)
## SINGLE LOCUS TREE-BASED
  ## NUCLEOTIDE DIVERSITY
    # start time of analysis from output of forward simulation
    tot_time=time.time()
    # start monitoring memory
    tracemalloc.start()
    # start time of calculation
    start_time = time.time()
    # calculate nucleotide diversity on individuals alive in final generation (whole population sample) using tree-based calculation
    mut_ts.diversity(sample_sets=pyslim.individuals_alive_at(mut_ts, 0))
    # end time of analysis
    end_time=time.time()
    # record the memory of the calculation IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    # stopping the library and adding recorded values to vectors
    tracemalloc.stop()
    total_time.append(end_time-tot_time)
    calc_time.append(end_time-start_time)
    fts_stat_type.append("Tree-based")
    fts_stat.append("div")
    fts_model.append("single_locus")
    dt.append("Tree Sequence")
    
  ## TAJIMAS D
    # start monitoring memory
    tracemalloc.start()
    tot_time=time.time()
    start_time = time.time()
    # calculate tajimas on individuals alive in final generation (whole population sample) using tree-based calculation
    mut_ts.Tajimas_D(sample_sets=pyslim.individuals_alive_at(mut_ts, 0))
    # record resource usage
    end_time=time.time()
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    tracemalloc.stop()
    total_time.append(end_time-tot_time)
    calc_time.append(end_time-start_time)
    fts_stat_type.append("Tree-based")
    fts_stat.append("tajimas d")
    fts_model.append("single_locus")
    dt.append("Tree Sequence")

## USING ALLELE-BASED CALCULATIONS ON TREE SEQUENCE OUTPUT
   ## NUCLEOTIDE DIVERSITY
    # start time of analysis from output of forward simulation
    tot_time = time.time()
    # sample individuals alive in the final generation
    samp_ts = mut_ts.simplify(samples=pyslim.individuals_alive_at(mut_ts, 0))
    # generate genotype matrix of those individuals
    gm=samp_ts.genotype_matrix()
    # convert genotype matrix to haplotyoe array 
    h= allel.HaplotypeArray(gm)
    # allele count for scikit.allel stats
    ac = h.count_alleles() 
    # ordered vector of mutation positions for scikit.allel stats
    mut_positions = [int(var.position+1) for var in samp_ts.variants()]
    # start monitoring memory
    tracemalloc.start()
    # start time for calculation
    start_time = time.time()
    # calculate nucleotide diversity using allele-based calculation
    allel.sequence_diversity(mut_positions, ac)
    # end time of analysis
    end_time=time.time()
    # record memory usage
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    # stopping the library and recording time, statistic and data type, and selection model
    tracemalloc.stop()
    total_time.append(end_time-tot_time)
    calc_time.append(end_time-start_time)
    fts_stat_type.append("Allele-based")
    fts_stat.append("div")
    fts_model.append("single_locus")
    dt.append("Tree Sequence")

   ## TAJIMAS D
    # start time of analysis from output of forward simulation
    tot_time = time.time()
    # sample individuals alive in the final generation
    samp_ts = mut_ts.simplify(samples=pyslim.individuals_alive_at(mut_ts, 0))
    # generate genotype matrix of those individuals
    gm=samp_ts.genotype_matrix()
    # convert genotype matrix to haplotyoe array
    h= allel.HaplotypeArray(gm)
    # allele count for scikit.allel stats
    ac = h.count_alleles()
    # ordered vector of mutation positions for scikit.allel stats
    mut_positions = [int(var.position+1) for var in samp_ts.variants()]
    # start monitoring memory
    tracemalloc.start()
    # start time for calculation
    start_time = time.time()
    # calculate tajimas D using allele-based calculation
    allel.tajima_d( ac, mut_positions)
    # end time of analysis
    end_time=time.time()
    # recording time, memory usage, statistic and data type, and selection model
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   # stopping the library
    tracemalloc.stop()
    calc_time.append(end_time-start_time)
    total_time.append(end_time-tot_time)
    fts_stat_type.append("Allele-based")
    fts_stat.append("tajimas d")
    fts_model.append("single_locus")
    dt.append("Tree Sequence")

## ALLELE-BASED CALCULATIONS ON CLASSICAL SIMULATION OUTPUT
   ## NUCLEOTIDE DIVERSITY
    # start time of analysis from output of forward simulation
    tot_time = time.time()
    # extract data for individuals in final generation and move to temporary file
    cmd="sed -n -e '/#OUT: 140000 140000 A/, /Individuals:/ p' {0}forward_single_mut{1}.txt > {0}final_gen.txt".format(tmpdir,i)
    os.system(cmd)
    # read in temporary file
    sts=pd.read_csv("{0}final_gen.txt".format(tmpdir), sep=" ", header=1, skiprows=4, skipfooter=1, names=["No", "ID", "mutType", "Pos", "s", "d", "pop", "Tick", "Count"])
    # order values by mutation position
    sts=sts.sort_values('Pos')
    # create alternate allele count by substracting current allele count from the number of total alleles (2 times the population size)
    sts['Alt']=(2*s_pop)-sts['Count']
    # create list of mutation positions for scikit allel function
    positions=list(sts["Pos"])
    # create allele count object from scikit allel function
    count=sts[["Count", "Alt"]].to_numpy()
    # start monitoring memory
    tracemalloc.start()
    #start time of calculation
    start_time = time.time()
    # calculate nucleotide diversity with allele-based calculation
    allel.sequence_diversity(positions, count)
    # end time of analysis
    end_time=time.time()
    # record resource usage and associated simualtion information
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    tracemalloc.stop()
    calc_time.append(end_time-start_time)
    total_time.append(end_time-tot_time)
    fts_stat_type.append("Allele-based")
    fts_stat.append("div")
    fts_model.append("single_locus")
    dt.append("Classical")
    
  ## TAJIMAS D
    # start time of analysis from output of forward simulation
    tot_time = time.time()
    # extract data for individuals in final generation and move to temporary file
    cmd="sed -n -e '/#OUT: 140000 140000 A/, /Individuals:/ p' {0}forward_single_mut{1}.txt > {0}final_gen.txt".format(tmpdir,i)
    os.system(cmd)
    # read in temporary file
    sts=pd.read_csv("{0}final_gen.txt".format(tmpdir), sep=" ", header=1, skiprows=4, skipfooter=1, names=["No", "ID", "mutType", "Pos", "s", "d", "pop", "Tick", "Count"])
    # order values by mutation position
    sts=sts.sort_values('Pos')
    # create alternate allele count
    sts['Alt']=20000-sts['Count']
    # create list of mutation positions for scikit allel function
    positions=list(sts["Pos"])
    # create allele count object from scikit allel function
    count=sts[["Count", "Alt"]].to_numpy()
    # start monitoring memory
    tracemalloc.start()
    # calculation start time
    start_time = time.time()
    # calculate tajimas d with allele-based calculation
    allel.tajima_d( count, positions)
    # end time of analysis
    end_time=time.time()
    # record resource usage and associated simualtion information
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    tracemalloc.stop()
    calc_time.append(end_time-start_time)
    total_time.append(end_time-tot_time)
    fts_stat_type.append("Allele-based")
    fts_stat.append("tajimas d")
    fts_model.append("single_locus")
    dt.append("Classical")
    
    
## ANALYSIS FOR MULTILOCUS MODEL (for code annotation see above)
for i in range(replicates):
    # load data
    slim_ts=tskit.load("{0}/treeseq_multi_{1}.trees".format(tmpdir,i)).simplify()
    
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)
## TREE-BASED CALCULATION ON TREE SEQUENCE
   ## NUCLEOTIDE DIVERSITY
    tracemalloc.start()
    
    tot_time=time.time()
    start_time = time.time()
    
    mut_ts.diversity(sample_sets=pyslim.individuals_alive_at(mut_ts, 0))
    
    end_time=time.time()
    
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    tracemalloc.stop()
    total_time.append(end_time-tot_time)
    calc_time.append(end_time-start_time)
    fts_stat_type.append("Tree-based")
    fts_stat.append("div")
    fts_model.append("multilocus")
    dt.append("Tree Sequence")

   ## TAJIMAS D
    tracemalloc.start()
    tot_time=time.time()
    start_time = time.time()
    
    mut_ts.Tajimas_D(sample_sets=pyslim.individuals_alive_at(mut_ts, 0))
  
    end_time=time.time()
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   
    tracemalloc.stop()
    total_time.append(end_time-tot_time)
    calc_time.append(end_time-start_time)
    fts_stat_type.append("Tree-based")
    fts_stat.append("tajimas d")
    fts_model.append("multilocus")
    dt.append("Tree Sequence")

    
## ALLELE-BASED CALCULATION ON TREE SEQUENCE
   ## NUCLEOTIDE DIVERSITY
    tot_time = time.time()
    
    samp_ts = mut_ts.simplify(samples=pyslim.individuals_alive_at(mut_ts, 0))
    gm=samp_ts.genotype_matrix()
    h= allel.HaplotypeArray(gm)
    
    ac = h.count_alleles()
    mut_positions = [int(var.position+1) for var in samp_ts.variants()]
    
    tracemalloc.start()
    start_time = time.time()
    allel.sequence_diversity(mut_positions, ac)
   
    end_time=time.time()
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    tracemalloc.stop()
    total_time.append(end_time-tot_time)
    calc_time.append(end_time-start_time)
    fts_stat_type.append("Allele-based")
    fts_stat.append("div")
    fts_model.append("multilocus")
    dt.append("Tree Sequence")

   
   ## TAJIMAS D
    tot_time = time.time()

    samp_ts = mut_ts.simplify(samples=pyslim.individuals_alive_at(mut_ts, 0))
    gm=samp_ts.genotype_matrix()
    h= allel.HaplotypeArray(gm)
    
    ac = h.count_alleles()
    mut_positions = [int(var.position+1) for var in samp_ts.variants()]
    
    tracemalloc.start()
    start_time = time.time()
    
    allel.tajima_d( ac, mut_positions)
  
    end_time=time.time()
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   
    tracemalloc.stop()
    calc_time.append(end_time-start_time)
    total_time.append(end_time-tot_time)
    fts_stat_type.append("Allele-based")
    fts_stat.append("tajimas d")
    fts_model.append("multilocus")
    dt.append("Tree Sequence")

## ALLELE-BASED CALCULATION ON CLASSICAL OUTPUT
   ## NUCLEOTIDE DIVERSITY
    tot_time = time.time()
    cmd="sed -n -e '/#OUT: 140000 140000 A/, /Individuals:/ p' {0}forward_multi_mut{1}.txt > {0}final_gen.txt".format(tmpdir,i)
    os.system(cmd)
    
    sts=pd.read_csv("{0}final_gen.txt".format(tmpdir), sep=" ", header=1, skiprows=4, skipfooter=1, names=["No", "ID", "mutType", "Pos", "s", "d", "pop", "Tick", "Count"])
    
    sts=sts.sort_values('Pos')
    sts['Alt']=20000-sts['Count']
    
    positions=list(sts["Pos"])
    count=sts[["Count", "Alt"]].to_numpy()
    
    tracemalloc.start()
    start_time = time.time()
    
    allel.sequence_diversity(positions, count)
    
    end_time=time.time()
    fts_memory.append(tracemalloc.get_traced_memory()[1])

    tracemalloc.stop()
    calc_time.append(end_time-start_time)
    total_time.append(end_time-tot_time)
    fts_stat_type.append("Allele-based")
    fts_stat.append("div")
    fts_model.append("multilocus")
    dt.append("Classical")
   
   ## TAJIMAS D
    tot_time = time.time()
    cmd="sed -n -e '/#OUT: 140000 140000 A/, /Individuals:/ p' {0}forward_multi_mut{1}.txt > {0}final_gen.txt".format(tmpdir,i)
    os.system(cmd)
    
    sts=pd.read_csv("{0}final_gen.txt".format(tmpdir), sep=" ", header=1, skiprows=4, skipfooter=1, names=["No", "ID", "mutType", "Pos", "s", "d", "pop", "Tick", "Count"])
    
    sts=sts.sort_values('Pos')
    sts['Alt']=20000-sts['Count']
    
    positions=list(sts["Pos"])
    count=sts[["Count", "Alt"]].to_numpy()
    
    tracemalloc.start()
    start_time = time.time()
    
    allel.tajima_d(count, positions)
    
    end_time=time.time()
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   
    tracemalloc.stop()
    calc_time.append(end_time-start_time)
    total_time.append(end_time-tot_time)
    fts_stat_type.append("Allele-based")
    fts_stat.append("tajimas d")
    fts_model.append("multilocus")
    dt.append("Classical")
  
# create table of data for analysis of all simulations
smdata = {'stat_type': fts_stat_type, 'data_type':dt, 'model': fts_model,'calc_time': calc_time,'total_time': total_time, 'stat': fts_stat, 'memory': fts_memory}
smbench_data = pd.DataFrame(smdata)
# export as txt file
smbench_data.to_csv("{0}/analysis_benchmarking.txt".format(tmpdir), index=False, header = True, sep = "\t")



    
    
    
    