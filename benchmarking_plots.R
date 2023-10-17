## Johnson et al 2023 - Population genetic simulation: Benchmarking frameworks for non-standard models of natural selection ##

library(data.table)
library(ggplot2)
library(lubridate)
library(ggpubr)
library(writexl)
library(scales)

setwd("~/2023_popgen_benchmarking/")

### BURN-IN 
  # read in msprime burn-in data
msprime_ts=fread("~/2023_popgen_benchmarking/burnin_benchmarking_msprime.txt")
  # add additional information regarding simulator, data type and final generation
msprime_ts[, simulator:="msprime"]
msprime_ts[, data_type:=ifelse(sim_type=="msprime_ts", "Tree Sequence","Tree Sequence\nwith Neutral Mutations")]
msprime_ts[, final_gen:=NA]

  # read in SLiM tree sequence (10N) data
slim_ts=fread("~/2023_popgen_benchmarking/burnin_benchmarking_slim_ts.txt")
  # format start and end time, calculate time for simulation
slim_ts[, start_time:=period_to_seconds(hms(start_time))]
slim_ts[,end_time:=period_to_seconds(hms(end_time))]
slim_ts[, time:=ifelse(end_time>start_time, end_time-start_time,((86400-start_time)+end_time)) ] # correct for simulations the ran past 12 am
  # add additional information regarding simulator, data type and final generation
slim_ts[, sim_type:="slim_ts"]
slim_ts[, simulator:="SLiM"]
slim_ts[, data_type:="Tree Sequence\n(10Ne)"]
slim_ts[, final_gen:=100000]

  # read in SLiM tree sequenc (checkCoalescence = T) data
slim_tscc=fread("~/2023_popgen_benchmarking/burnin_benchmarking_slim_tscc.txt")
  # format start and end time, calculate time for simulation
slim_tscc[, start_time:=period_to_seconds(hms(start_time))]
slim_tscc[,end_time:=period_to_seconds(hms(end_time))]
slim_tscc[, time:=ifelse(end_time>start_time, end_time-start_time,((86400-start_time)+end_time)) ]
  # add additional information regarding simulator and data type
slim_tscc[, sim_type:="slim_tscc"]
slim_tscc[, simulator:="SLiM"]
slim_tscc[, data_type:="Tree Sequence\n(Coalesced)"]

  # read in classical SLiM burn-in data
slim_mut=fread("~/2023_popgen_benchmarking/burnin_benchmarking_slim_mut.txt")
  # format start and end time, calculate time for simulation
slim_mut[, start_time:=period_to_seconds(hms(start_time))]
slim_mut[,end_time:=period_to_seconds(hms(end_time))]
slim_mut[, time:=ifelse(end_time>start_time, end_time-start_time,((86400-start_time)+end_time))]
  # add additional information regarding simulator, data type and final generation
slim_mut[, sim_type:="slim_mut"]
slim_mut[, simulator:="SLiM"]
slim_mut[, data_type:="Classical\n(10Ne)"]
slim_mut[, final_gen:=100000]

# compile into single data object
data=rbind(msprime_ts,slim_ts[,.(time, memory, diversity,sim_type, simulator, data_type, final_gen)],slim_tscc[,.(time, memory, diversity,sim_type, simulator, data_type, final_gen)], slim_mut[,.(time, memory, diversity, sim_type, simulator, data_type,final_gen)])
data[, memory:=memory/1e6]  ## convert memory to MB

# calculate mean, variance, min and max of each resource for each simualtion type
data[, `:=` (mean_time = mean(time), mean_memory=mean(memory), var_time=var(time), var_mem=var(memory), max_time=max(time), min_time=min(time), max_mem=max(memory), min_mem=min(memory), mean_div=mean(diversity), max_div=max(diversity), min_div=min(diversity),var_div=var(diversity), var_fg=var(final_gen), mean_fg=mean(final_gen), min_fg=min(final_gen), max_fg=max(final_gen)), by="sim_type"]

# table of just summary values
averages = unique(data[,.(simulator, data_type, mean_time, mean_memory, var_time, var_mem, max_mem, min_mem, max_time, min_time, mean_div, min_div, max_div, var_div, mean_fg, var_fg, min_fg, max_fg)])

## Calculate significance usign t-tests
  ## significant difference in memory  between msprime tree sequence and tree sequence with neutral mutation
t.test(data[sim_type=="msprime_ts", memory], data[sim_type=="msprime_mut", memory])
  ## significcant difference in time  between msprime tree sequence and tree sequence with neutral mutation
t.test(data[sim_type=="msprime_ts", time], data[sim_type=="msprime_mut", time])

  ## Is diversity significantly different from neutral expectation
theta=4*10000*1e-7
expected_div=theta/(1+theta)
t.test(data[sim_type=="msprime_ts", diversity], mu=expected_div) # msprime coalescent with just tree sequence
t.test(data[sim_type=="msprime_mut", diversity], mu=expected_div) # msprime coalescent with neutral mutations
t.test(data[sim_type=="slim_ts", diversity], mu=expected_div) # SLiM tree sequence (10Ne)
t.test(data[sim_type=="slim_tscc", diversity], mu=expected_div) # SLiM tree sequence (checkCoalescence = T)
t.test(data[sim_type=="slim_mut", diversity], mu=expected_div) # SLiM classical (10Ne)


## FIGURE 3
mem = ggplot(data, aes(x=data_type))+
  geom_boxplot(aes( y=memory, fill=sim_type))+
  labs(x="Simulation Type",y="Memory Usage (MB)") +
  theme_bw()+
  theme(legend.position = "none")+ 
  facet_wrap("simulator", scale="free", drop=TRUE)+
  scale_fill_manual(values = c("#C77CFF", "#F8766D", "#00BFC4", "#F8766D", "#F8766D"))

time = ggplot(data,aes(x=data_type))+
  geom_boxplot(aes(y=time, fill=sim_type))+
  labs(x="",y="Time (seconds)")+
  theme_bw() +
  theme(legend.position = "none")+ 
  facet_wrap("simulator", scale="free", drop=TRUE)+ 
  scale_fill_manual(values = c("#C77CFF", "#F8766D", "#00BFC4", "#F8766D", "#F8766D"))
ggexport(burnin, filename="figure_3.pdf", width=7.5, height=7)
ggsave(plot=burnin, filename="figure_3.jpg", width=7.5, height=7)


## FIGURE 4
div=ggplot(data[sim_type!="msprime_ts"], aes(x=data_type)) + 
  geom_boxplot(aes(y = diversity/expected_div, fill=sim_type))+
  labs(x="Simulation Type",y="Relative Diversity") +
  theme_bw()+theme(legend.position = "none")+ 
  facet_grid(~simulator, drop=TRUE, scale='free_x', space='free_x') +
  scale_fill_manual(values = c("#C77CFF", "#F8766D", "#00BFC4", "#F8766D", "#F8766D")) 
ggexport(div, filename="figure_4.pdf", height = 4, width=8)
ggsave(plot=div, filename="figure_4.jpg", height = 4, width=8)

### FORWARD SIMUALTION OF SELECTION

  # read in resource usage data
for_t_ts=fread("~/2023_popgen_benchmarking/forward_benchmarking_time_ts.txt")
for_m_m=fread("~/2023_popgen_benchmarking/forward_mut_benchmarking_time_multi.txt")
for_m_s=fread("~/2023_popgen_benchmarking/forward_mut_benchmarking_time_single.txt")
forward_time=rbind(for_t_ts,for_m_m,for_m_s)

forward_mem=fread("~/2023_popgen_benchmarking/forward_benchmarking_mem.txt")
   
  #calculate mean, variance, minimum and maximum of resurces for each simulation
forward_time[, `:=` (mean_time = mean(time), var_time=var(time),  max_time=max(time), min_time=min(time)), by=c("sim_type", "model")]
forward_mem[, `:=` ( mean_memory=mean(memory),  var_mem=var(memory), max_mem=max(memory), min_mem=min(memory)), by=c("sim_type", "model")]
  # add additional labels
forward_time[, sim_type:=ifelse(sim_type=="mut", "Classical", "Tree Sequence"), by="sim_type"]
forward_time[, model:=ifelse(model=="single_locus", "Single Locus", "Multilocus"), by="model"]
forward_mem[, sim_type:=ifelse(sim_type=="mut", "Classical", "Tree Sequence"), by="sim_type"]
forward_mem[, model:=ifelse(model=="single_locus", "Single Locus", "Multilocus"), by="model"]
  # formulate summary table
forward=merge(unique(forward_mem[, .(sim_type, model, mean_memory, var_mem, max_mem, min_mem)]), unique(forward_time[, .(sim_type, model, mean_time, var_time, max_time, min_time)]), by=c("sim_type", "model"))

## FIGURE 5
plot1 = ggplot(forward_mem)+
  geom_boxplot(aes(x=sim_type, y=memory, fill=sim_type))+
  facet_wrap(~factor(model, levels=c("Single Locus", "Multilocus")))+theme_bw()+ 
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="Memory (MB)", x="Simulation Type")+
  theme(legend.position = "none")+
  scale_y_log10()

plot2 = ggplot(forward_time)+
  geom_boxplot(aes(x=sim_type, y=time, fill=sim_type))+
  theme_bw()+
  facet_wrap(~factor(model, levels=c("Single Locus", "Multilocus")))+ 
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="Time (seconds)", x="Simulation Type")+
  scale_y_log10()+
  theme(legend.position = "none", axis.title.x.bottom = element_blank())

forwardp= ggarrange(plot2, plot1, ncol = 1, labels = c("A", "B"))
ggexport(forwardp, filename="figure_5.pdf")
ggsave(plot=forwardp, filename="figure_5.jpg", height = 7, width=7.5)

### ANALYSIS 
  # read in resource usage data
analysis=fread("~/2023_popgen_benchmarking/analysis_benchmarking.txt")
  # convert memory to MB
analysis[, memory:=as.numeric(memory)/1e6] 
  # add additional labels
analysis[, model:=ifelse(model=="single_locus", "Single Locus", "Multilocus"), by="model"]
analysis[, stat:=ifelse(stat=="div", "Nucleotide Diversity", "Tajima's D"), by="stat"]
analysis[, x_lab:=paste(data_type, "\n", stat_type), by=c("data_type", "stat_type")]
  # calculate mean, variance, minimum and maximum for each statistics, simulation and calculation type
analysis[, `:=` (mean_calc_time = mean(calc_time),mean_total_time = mean(total_time), mean_memory=mean(memory), var_total_time=var(total_time),var_calc_time=var(calc_time), var_mem=var(memory), max_total_time=max(total_time), min_total_time=min(total_time),max_calc_time=max(calc_time), min_calc_time=min(calc_time), max_mem=max(memory), min_mem=min(memory)), by=c("stat_type","data_type", "model", "stat")]
  # summary values for analysis results
analysis_sum = unique(analysis[, .(stat_type, data_type, model, stat, mean_memory, var_mem, max_mem, min_mem, mean_calc_time, var_calc_time, min_calc_time, max_calc_time, mean_total_time, var_total_time, min_total_time, max_total_time)])

## FIGURE 6
  # memory usage
plot3.1 = ggplot(analysis[stat=="Nucleotide Diversity"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_boxplot(aes(x=stat_type, y=memory, fill=stat_type))+ 
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  theme_bw()+coord_cartesian(y=c(0, 2.25))+
  facet_grid(model~stat)+
  theme( axis.title.y.left = element_blank(),legend.position = "none", axis.title.x.bottom = element_blank(), 
         strip.text.y = element_blank())
plot3.2 = ggplot(analysis[stat=="Tajima's D"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_boxplot(aes(x=stat_type, y=memory,fill=stat_type))+
  theme_bw()+ scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  facet_grid(model~stat)+
  coord_cartesian(y=c(0, 2.25))+
  theme(strip.text.y = element_blank(),legend.position = "none", axis.title.y.left = element_blank(), 
        axis.title.x.bottom = element_blank())
plot3=ggarrange(plot3.1, plot3.2, ncol=1)
plot3=annotate_figure(plot3, left = "Memory (MB)")

  # calculation time
plot4.1 = ggplot(analysis[stat=="Nucleotide Diversity"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_boxplot(aes(x=stat_type, y=calc_time,fill=stat_type))+
  theme_bw()+ 
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  facet_grid(model~stat)+
  coord_cartesian(y=c(0, 1.5))+
  theme( axis.title.y.left = element_blank(),legend.position = "none", axis.title.x.bottom = element_blank(), 
         strip.text.y = element_blank())
plot4.2 = ggplot(analysis[stat=="Tajima's D"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_boxplot(aes(x=stat_type, y=calc_time,fill=stat_type))+
  theme_bw()+ 
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  facet_grid(model~stat)+ 
  coord_cartesian(y=c(0, 1.5))+
  theme(strip.text.y = element_blank(),legend.position = "none", axis.title.y.left  = element_blank(), 
        axis.title.x.bottom = element_blank())
plot4=ggarrange(plot4.1, plot4.2, ncol=1)
plot4=annotate_figure(plot4, left = "Calculation time (seconds)")


  # total time
plot5.1 = ggplot(analysis[stat=="Nucleotide Diversity"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_boxplot(aes(x=stat_type, y=total_time,fill=stat_type))+
  theme_bw()+ scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  facet_grid(model~stat)+
  coord_cartesian(y=c(0, 1.5))+
  theme( axis.title.y.left = element_blank(),legend.position = "none", axis.title.x.bottom = element_blank())
plot5.2 = ggplot(analysis[stat=="Tajima's D"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_boxplot(aes(x=stat_type, y=total_time,fill=stat_type))+
  theme_bw()+ 
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  facet_grid(model~stat)+
  coord_cartesian(y=c(0, 1.5))+
  labs(y="Total time (seconds)", x="")+
  theme(legend.position = "none", axis.title.y.left  = element_blank(), axis.title.x.bottom = element_blank())
plot5=ggarrange(plot5.1, plot5.2, ncol=1)
plot5=annotate_figure(plot5, left="Total time (seconds)")

# compile to one plot
stat_an= ggarrange(plot3,plot4,plot5,  nrow = 1, labels = c("A", "B", "C"))
stat_an=annotate_figure(stat_an, bottom = "Calculation Type")
ggexport(stat_an, filename="figure_6.pdf", height=6, width=9)
ggsave(plot=stat_an, filename="figure_6.jpg", height = 6, width=9)

## export to xl
write_xlsx(list( burnin=averages, forward=forward, analysis=analysis_sum), path="benchmark_summary.xlsx",col_names=TRUE)


