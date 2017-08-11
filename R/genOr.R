
genOr <- function(threshold, kME, cutVal){
	if(missing(cutVal)){ cutVal = TRUE }
	premergeKME <- kME
	receiving <- data.frame(matrix(0, nrow=nrow(premergeKME), ncol=ncol(premergeKME)))

for(k in 1:nrow(premergeKME)) {
	shelf <- data.frame(-sort(premergeKME[k,]))
	s <- sum(shelf[1,]>threshold)
	if(s>=1){
		col <-colnames(shelf[1:s])
		receiving[k, 1:length(col)] <- col
}
}


#count number of modules each gene is assigned to
rownames(premergeKME) -> rownames(receiving)
Out <- receiving
count <- data.frame(matrix(0, (nrow(Out)), 1))

for(k in 1:nrow(Out)){
	count[k,] <- sum(Out[k,] != 0)
	}
rownames(count) <- rownames(Out)
colnames(count) <- c("Module Memberships")
count$threshold <- threshold
Mod <- cbind(count, Out)
tab <- table(Mod$"Module Memberships")

if(cutVal == TRUE) { Mod <- Mod[,1:(length(tab)+1)]}

list(GenOrOutput <- Mod, Tabled_Membership <- tab)
}
