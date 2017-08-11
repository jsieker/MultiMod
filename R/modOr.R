#`
#`

modOr <- function(kME, cap, threshold){
if(missing(cap)){cap <- 5}
if(missing(threshold)){threshold <- 0.75}
if(missing(kME)){stop()}
premergeKME <- kME
premade <- matrix(nrow = nrow(premergeKME), ncol = (2 * cap))
Top5 <- data.frame(premade)
colnames(Top5) <- c(1:(2*cap))
row.names(Top5) <- row.names(premergeKME)


for(k in 1:nrow(premergeKME)) {
	shelf <- data.frame(-sort(premergeKME[k,])[1:cap])
	tempnames <- colnames(shelf)
	combined <- rbind(shelf, tempnames)
	for(j in 1:cap){
	ifelse(combined[1,j]<threshold, combined[3,j]<-1, combined[3,j]<-0)
	ifelse(combined[3,j]==1, combined[,j]<- "NONE", combined[3,j]<-NA)
	}
	Top5[k, c(TRUE, FALSE)] <- combined[1,]
	Top5[k, !c(TRUE, FALSE)] <- combined[2,]
}


GeneModuleDF <- data.frame(Top5[,!c(TRUE, FALSE)])
ModuleNames <- sort(colnames(premergeKME))
GeneModuleList <- list()
for(p in 1:nrow(GeneModuleDF)){GeneModuleList[[p]]<-0}
names(GeneModuleList) <- row.names(GeneModuleDF)
OutOut <- list()

if(cap > 1){
for(i in 1:length(ModuleNames)) {		#this one is going to need to be redesigned for other caps
	as.character(ModuleNames[i]) -> y
	set <- apply(GeneModuleDF, 1, function(x) { any(x == y)})
	set2 <- set[set==TRUE]
	shelf <- names(set2)
	if(length(shelf)!=0){
	OutOut[[i]] <- shelf
	} else {OutOut[[i]] <- "<NA>"}
}
}

if(cap == 1){
	for(i in 1:length(ModuleNames)) {		#this one is going to need to be redesigned for other caps
	as.character(ModuleNames[i]) -> y
	set <- apply(GeneModuleDF, 1, function(x) { any(x == y)})
	set2 <- row.names(Top5)[set]
	shelf <- set2
	if(length(shelf)!=0){
	OutOut[[i]] <- shelf
	} else {OutOut[[i]] <- "<NA>"}
}
}




names(OutOut) <- ModuleNames
xx <- lapply(OutOut, unlist)
max <- max(sapply(xx, length))
exportable <- data.frame(t(do.call(rbind, lapply(xx, function(z)c(z, rep(NA, max-length(z)))))))

}
