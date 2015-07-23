################################################
# A script to run the rehh genome scan for EHH #
################################################
  
#Simon Renny-Byfield, UC Davis, July 2015
  
# load in the two command line arguments
# arg1 = fastPHASE haplotypes
# agr2 = map info file

args <- commandArgs(trailingOnly = TRUE)
print(args)

############################
# now perform the analysis #
 ################################
 
print(paste0("loading file:",args[1]))
# make the hap object from the fastPHASE input
hap<-data2haplohh(hap_file=args[i],map_file=args[2])
print("Performing Scan")
# scan genome wide
res<-scan_hh(hap)
save(res,file=paste0(args[1],"rehh.RData"))
