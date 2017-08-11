
modPlots <- function(expr, modOrOut, ord){

if(missing(ord) | ord == "intensity"){
  flag <- 0 #intensity-wise
  ord <- 0
  }
if(ord=="alpha") {
  flag <- 1
  }
if(length(ord)==ncol(exportable)) {
  flag <- 2
  }
if(flag != (0 | 1 | 2)) {
  flag <- 0
  }

exportable <- modOrOut
indent = 0
nPC <- 1
alignRecognizedValues = c("", "along average")
maxVarExplained = 10
nVarExplained = min(nPC, maxVarExplained)
modlevels = levels(factor(colnames(exportable)))
PrinComps = data.frame(matrix(NA, nrow = dim(expr)[[1]], ncol = length(modlevels)))
averExpr = data.frame(matrix(NA, nrow = dim(expr)[[1]], ncol = length(modlevels)))
varExpl = data.frame(matrix(NA, nrow = nVarExplained, ncol = length(modlevels)))
validMEs = rep(TRUE, length(modlevels))
validAEs = rep(FALSE, length(modlevels))
isPC = rep(TRUE, length(modlevels))
isHub = rep(FALSE, length(modlevels))
validColors = colnames(exportable)
names(PrinComps) = paste(moduleColor.getMEprefix(), modlevels, sep = "")
names(averExpr) = paste("AE", modlevels, sep = "")

for (i in c(1:length(modlevels))) {
  shelf <- exportable[,i]
  shelf <- shelf[which(!is.na(shelf))]
  if(!is.na(shelf[1])){
  datModule <- as.matrix(t(expr[, (colnames(expr) %in% shelf)]))
  n = dim(datModule)[1]
  p = dim(datModule)[2]
  pc = try({
    datModule = t(scale(t(datModule)))
    svd1 = svd(datModule, nu = min(n, p, nPC), nv = min(n,
                                                        p, nPC))
    veMat = cor(svd1$v[, c(1:min(n, p, nVarExplained))],
                t(datModule), use = "p")
    varExpl[c(1:min(n, p, nVarExplained)), i] = rowMeans(veMat^2,
                                                         na.rm = TRUE)
    svd1$v[, 1]
  }, silent = TRUE)

  PrinComps[, i] = pc
  ae = try({
    if (isPC[i])
      scaledExpr = scale(t(datModule))
    averExpr[, i] = rowMeans(scaledExpr, na.rm = TRUE)
    if (align == "along average") {
      if (verbose > 4)
        printFlush(paste(spaces, " .. aligning module eigengene with average expression."))
      corAve = cor(averExpr[, i], PrinComps[, i],
                   use = "p")
      if (!is.finite(corAve))
        corAve = 0
      if (corAve < 0)
        PrinComps[, i] = -PrinComps[, i]
    }
    0
  }, silent = TRUE)

  validAEs[i] = !(class(ae) == "try-error")
  } else { PrinComps[, i] <- 0}
}

row.names(PrinComps) <- row.names(expr)


plotOutput <- list()
for (k in 1:ncol(PrinComps)){
col <- k
row1 <- row.names(PrinComps)#x axis
row2 <- PrinComps[,col]
rrs <- data.frame(row1, row2)
if(flag==0) { rrs_s <- rrs[order(rrs$row2),]; ord <- rrs_s$row1 } #ME-ordered
if(flag==1) { rrs_s <- rrs; ord <- rrs$row1} #alphabetical
if(flag==2) { rrs_s <- rrs; } #user defined
if(sum(rrs$row2) != 0){
p <-ggplot(rrs_s, aes(row1, row2))
x <- p +geom_bar(stat = "identity", fill = "navy blue") +theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(paste("Modular Relative Expression", modlevels[col], sep = "--")) +
  ylab("Expression of genes in module") + xlab("Sample") + scale_x_discrete(limits = ord)
plotOutput[[k]] <- x
} else {plotOutput[[k]] <- 0}
}

plotOutput


} #end function

