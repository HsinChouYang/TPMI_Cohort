library("admixtools")
library("tidyverse")

load("J:/TPMI_F4scpre/Sample600_chr6.RData")

# Start with a graph with 0 admixture events, increase up to 3, and stop after 10 generations of no improvement
pops = dimnames(f2_blocks)[[1]]
re_ini=T
score=NULL
for(i in 1:1000){
	if(re_ini==T){
		initgraph = random_admixturegraph(pops, 0)
	}else{
		initgraph = winner$graph[[1]]
	}
	res = find_graphs(f2_blocks, initgraph = initgraph, stop_gen = 100, stop_gen2 = 15, max_admix = 3)
	res %>% slice_min(score)
	winner = res %>% slice_min(score, with_ties = FALSE)
	temp_score=winner$score[[1]]

	score=c(score,temp_score)
	if(temp_score==min(score)){
		re_ini=F
	}else{
		re_ini=T
	}
	print("==============================")
	print(paste(re_ini,"/",min(score)))
	print("==============================")
}

