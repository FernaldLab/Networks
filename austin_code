
#################################################################################
### load raw TPM data
#################################################################################
rm(list=ls());
options(stringsAsFactors=F);
source('/Volumes/fishstudies-1/_code/preProcATH-for_web_noVSN.R');

expressionDIR = '~/Documents/kallisto_expr_cutadaptJan27/'  
setwd(expressionDIR);

tpm0 = read.table(paste0(expressionDIR, 'TPM_raw.tsv'), header=T, row.names=1);
rg0 = read.table('~/Documents/rg_withFirst30Beh.tsv', header=T, row.names=1, sep='\t');

# make sure the samples match up before doing this
names(tpm0) = rownames(rg0);

#################################################################################
### define which samples to use
#################################################################################
samplesAll = 1:ncol(tpm0);
# all except ND
samplesExceptND = grep('^ND', names(tpm0), invert=T);
# only females
samplesFemales = grep('^FEM', names(tpm0));
# only males
samplesMales = grep('^FEM', names(tpm0), invert=T);
# only females and D males
samplesFemalesDmales = grep('^FEM|^D', names(tpm0));
# only ascenders and doms
samplesASCandD = grep('^A|^D', names(tpm0));
# only D and NDs
samplesOnlyD_ND = grep('^A|^FEM', names(tpm0), invert=T);
# only ASC and NDs
samplesOnlyASC_ND = grep('^A|^ND', names(tpm0));
# 
samplesD = grep('^D', names(tpm0));

SUBSET = samplesMales;
TPM = tpm0[, SUBSET];
INFO0 = cbind(Sample=rownames(rg0), rg0);
INFO0 = INFO0[SUBSET, ];


# ############################################################################################################
# ### remove genes with zero variance, too many zeros, or low average TPM



#toRemove = c(
#             'D_B6_lib1-ACAGTG-Mar2015-2',
#             'FEM_B6_lib3-ACAGTG-Dec2015-6'
 #            ,'ND_B6_lib3-CAGATC-Dec2015-1'
 #             ,'D_B6_lib4-GCCAAT-Dec2015-5'
#              ,'ND_B5_lib1-CCGTCC-Jan2015-2'
#             ,'D_B6_lib3-CGATGT-Dec2015-4'
#              ,'D_B5_lib2-TGACCA-Mar2015-1'
#              ,'FEM_B7_lib1-GTGAAA-Jan2015-2'
#);

INFO0 = cbind(Sample=rownames(rg0), rg0);
INFO0 = INFO0[SUBSET, ];

## pretty good for all males, followed by 0 removal, sd2:
toRemove = c(
  'ND_B6_lib3-CAGATC-Dec2015-1'
 # ,'D_B6_lib4-GCCAAT-Dec2015-5'
 # ,'ND_B5_lib4-GTGAAA-Dec2015-3'
 # ,'D_B6_lib4-CTTGTA-Dec2015-6'
 # ,'ND_B6_lib3-CTTGTA-Dec2015-3'
 # ,'D_B6_lib3-CGATGT-Dec2015-4'
  #, 'D_B7_lib4-GCCAAT-Dec2015-3'
  #, 'D_B6_lib4-CGATGT-Dec2015-7'
      #        
)
## pretty good for DvsND: remove ND_B6_lib3-CAGATC-Dec2015-1, D_B6_lib4-GCCAAT-Dec2015-5, D_B6_lib4-CTTGTA-Dec2015-6

 


DAT0 = TPM[, -match(toRemove, names(TPM))]
INFO0 = INFO0[-match(toRemove, rownames(INFO0)), ];
# DAT0 = TPM
# remove zero variance genes
DAT0 = DAT0[apply(DAT0,1,var) > 0, ];



DAT0 = log2(DAT0 + 1);



# .removeHighTPMhelper = function (tpm, numSd=1) {
#   highExprThresh = min(apply(tpm, 2, max, na.rm=T)) - 1;
#   tpm[which(tpm > highExprThresh, arr.ind=T)] = NA;
#   newMaxes = apply(tpm, 2, max, na.rm=T);
#   maxsd = sd(newMaxes) * numSd;
#   checkMax = max(newMaxes) > (mean(newMaxes) + maxsd);
#   return(list(test=checkMax, new=tpm, maxes=newMaxes));
# }
# 
# .removeHighTPM = function (tpm, numSd=1) {
#   out = .removeHighTPMhelper(tpm, numSd=numSd);
#   run = 1;
#   while (out$test) {
#     print(run)
#     out = .removeHighTPMhelper(out$new, numSd=numSd);
#     run = run+1;
#   } 
#   return(out);
# }
# 
# DAT0 = .removeHighTPM(DAT0,2)$new;


# set limit for number of zeros allowed
#zlim = floor(ncol(DAT0)/3);
zlim = floor(min(table(INFO0$Condition)/3));
#zlim = 0;
zs = .getZerosAcrossConditions(DAT0, INFO0$Condition); apply(zs, 2, table);

rmZeros = which((zs$zCounts > zlim) | (zs$zCounts > 1 & !zs$sameCond));
#rmZeros = which(zs$zCounts > zlim);

DAT0z = DAT0[-rmZeros, ];
summary(t(DAT0[rownames(DAT0)=='ptgfr',]));
summary(t(DAT0z[rownames(DAT0z)=='ptgfr',]));

dev.off(); 
par(mfrow=c(1,2))
xcv = apply(DAT0z, 1, cv); hist(xcv); summary(xcv)
cvcut = quantile(xcv, .01)
DAT0z = DAT0z[xcv > cvcut, ];
hist(apply(DAT0z, 1, cv)); summary(apply(DAT0z, 1, cv))


# ptgfr       
# Min.   :0.0000  
# 1st Qu.:0.8684  
# Median :1.1356  
# Mean   :1.1633  
# 3rd Qu.:1.4389  
# Max.   :2.2468 

# mtpm = apply(DAT0z, 1, mean);
# centile = 0.05;
# DAT0z = DAT0z[mtpm >= quantile(mtpm,centile), ];

#DAT0z[which(DAT0z < 1, arr.ind=T)] = NA;

# medtpm = apply(DAT0z, 1, median, na.rm=T); summary(medtpm)
# #medCut = median(medtpm);
# medCut = 1
# #medCut = 4.8
# DAT0z = DAT0z[(medtpm > medCut) & !is.na(medtpm), ];

#DAT0z = DAT0
DATout = preProcess(datIN=DAT0z, 
                    removeOutlierProbes=T, deviate=2.5, 
                    removeTooManyNAs=T, probe_thresh=NULL, 
                    sample_thresh=NULL, removeOutlierSamples=T, IACthresh=2, 
                    Qnorm=T);


dev.off()
par(mfrow=c(3,3));
FACTNAMES = c('Condition','Tank','Lib.constr.date','RNAseq.date');
.computeAndPlotFactorEffectsOnMeanExpr(DAT0z, INFO0, 1:ncol(DAT0z), match(FACTNAMES, names(INFO0)), '');
FACTNAMES = c('Condition','Tank','LibSeq');
.computeAndPlotFactorEffectsOnMeanExpr(DAT0z, INFO0, 1:ncol(DAT0z), match(FACTNAMES, names(INFO0)), '');
FACTNAMES = c('Condition','LibSeqTank');
.computeAndPlotFactorEffectsOnMeanExpr(DAT0z, INFO0, 1:ncol(DAT0z), match(FACTNAMES, names(INFO0)), '');

preNormPostOutliersDAT = DATout$data_removedOutlierSamples;
preNormPostOutliersINFO = INFO0[match(names(DATout$data_removedOutlierSamples), rownames(INFO0)), ];
FACTNAMES = c('Condition','Tank','Lib.constr.date','RNAseq.date');
.computeAndPlotFactorEffectsOnMeanExpr(preNormPostOutliersDAT, preNormPostOutliersINFO, 1:ncol(preNormPostOutliersDAT), match(FACTNAMES, names(INFO0)), '');
FACTNAMES = c('Condition','Tank','LibSeq');
.computeAndPlotFactorEffectsOnMeanExpr(preNormPostOutliersDAT, preNormPostOutliersINFO, 1:ncol(preNormPostOutliersDAT), match(FACTNAMES, names(INFO0)), '');
FACTNAMES = c('Condition','LibSeqTank');
.computeAndPlotFactorEffectsOnMeanExpr(preNormPostOutliersDAT, preNormPostOutliersINFO, 1:ncol(preNormPostOutliersDAT), match(FACTNAMES, names(INFO0)), '');

sort(table(paste(preNormPostOutliersINFO$Condition, preNormPostOutliersINFO$Lib.constr.date)))
sort(table(paste(preNormPostOutliersINFO$Condition, preNormPostOutliersINFO$RNAseq.date)))
sort(table(paste(preNormPostOutliersINFO$Condition, preNormPostOutliersINFO$Tank)))
sort(table(paste(preNormPostOutliersINFO$Condition, preNormPostOutliersINFO$LibSeq)))
sort(table(paste(preNormPostOutliersINFO$Condition, preNormPostOutliersINFO$LibSeqTank)))


DAT = DATout$data_Qnorm; 
INFO = preNormPostOutliersINFO;
FACTS = match(c('Condition','Tank','Lib.constr.date','RNAseq.date'), names(INFO));
#FACTS = match(c('Condition','Tank','LibSeq'), names(INFO));
#FACTS = match(c('Condition','LibSeqTank'), names(INFO));

BATCH = 'Lib.constr.date';
#BATCH = 'LibSeq';
BATCH = match(BATCH, names(INFO));
SAMPLE = match('Sample', names(INFO));
COVAR = match(c('Condition'), names(INFO));
x = .runComBat(DAT,INFO,SAMPLE,BATCH,COVAR,impute=T,prior.plots=F);
#dev.off();
.computeAndPlotFactorEffectsOnMeanExpr(x, INFO, 1:ncol(x), FACTS, '');
#xos = outlierSamples(x);

BATCH = 'Tank';
BATCH = match(BATCH, names(INFO));
SAMPLE = match('Sample', names(INFO));
COVAR = match('Condition', names(INFO));
xx = .runComBat(x,INFO,SAMPLE,BATCH,COVAR,impute=F,prior.plots=F);
#dev.off();
.computeAndPlotFactorEffectsOnMeanExpr(xx, INFO, 1:ncol(xx), FACTS, '');
#xxos = outlierSamples(xx);

BATCH = 'RNAseq.date';
BATCH = match(BATCH, names(INFO));
SAMPLE = match('Sample', names(INFO));
COVAR = match('Condition', names(INFO));
xxx = .runComBat(xx,INFO,SAMPLE,BATCH,COVAR,impute=F,prior.plots=F);
#dev.off();
.computeAndPlotFactorEffectsOnMeanExpr(xxx, INFO, 1:ncol(xxx), FACTS, '');
#xxxos = outlierSamples(xxx);


rownames(xxx) = rownames(DAT)
############################################################################################################
### functions 
############################################################################################################

.computeAndPlotFactorEffectsOnExpr = function (TPM, INFO0, SUBSET, INFOcols, MAIN='', ...) {
  DAT = TPM[, SUBSET];
  INFO = INFO0[SUBSET, ];
  DAT = DAT[apply(DAT,1,var) > 0, ];
  expr = c();
  for (j in 1:ncol(DAT)) { expr = c(expr, DAT[, j]) }
  kwpvals = c();
  for (f in INFOcols) {
    cat(paste(' ',names(INFO)[f],sep=''));
    fact = c();
    for (j in 1:ncol(DAT)) { fact = c(fact, rep(INFO[j, f], nrow(DAT))) };
    kwpvals = c(kwpvals, kruskal.test(as.numeric(expr), as.factor(fact))$p.value);
  }
  names(kwpvals) = names(INFO)[INFOcols];
  barplot(-log10(kwpvals), ylab='-log10(KW.pval)',main=paste(MAIN,' (n=',ncol(DAT),')',sep=''), ...);
  abline(h=c(-log10(.05), -log10(.01)), col='red');
}

.computeAndPlotFactorEffectsOnMeanExpr = function (TPM, INFO, SUBSET, INFOcols, MAIN='', ...) {
  DAT = TPM[, SUBSET];
  INFO = INFO[SUBSET, ];
  DAT = DAT[apply(DAT,1,var) > 0, ];
  meanexpr = apply(DAT, 2, mean, na.rm=T);
  factors = paste0(paste0('factor(',paste0('INFO[,',INFOcols,']' ),')' ), collapse='+');
  lmformula = paste0('meanexpr~', factors);print(lmformula)
  pvals = as.vector(na.omit(anova(lm(as.formula(lmformula)))$"Pr(>F)"));
  print(anova(lm(as.formula(lmformula))));
  names(pvals) = names(INFO)[INFOcols];
  barplot(-log10(pvals), ylab='-log10(pval)',main=paste(MAIN,' (n=',ncol(DAT),')',sep=''), ...);
  abline(h=c(-log10(.05), -log10(.01)), col='red');
}

source('/Volumes/fishstudies/_code/ComBat.R');
library(impute);

# rows of INFO MUST must match cols of DAT
.runComBat = function (DAT, INFO, SAMPLEcol, BATCHcol, COVcol=NULL, impute=F, ...) {
  if (impute) { DAT = as.data.frame(impute.knn(as.matrix(DAT))$data) }
  # format expression data and write to file
  datToWrite = as.data.frame(DAT);
  datToWrite = as.data.frame(cbind(gene=rownames(datToWrite), datToWrite));
  write.table(datToWrite, file='datForComBat.txt', quote=F, sep='\t', row.names=F);
  
  if (is.null(COVcol)) {
    inds = c(SAMPLEcol, BATCHcol);
  } else if (is.numeric(COVcol)) {
    inds = c(SAMPLEcol, BATCHcol, COVcol);
  } else {
    stop('arg COVcol must be NULL or numeric');
  }
  print(inds);
  infoToWrite = INFO[, inds];
  infoToWrite = as.data.frame(cbind(names(DAT), infoToWrite));
  names(infoToWrite)[1:3] = c('Array name', 'Sample name', 'Batch');
  if (!is.null(COVcol) & is.numeric(COVcol)) {
    names(infoToWrite)[4:ncol(infoToWrite)] = paste('Covariate ', 1:length(COVcol), sep='');	
  }
  write.table(infoToWrite, file='infoForComBat.txt', quote=F, sep='\t', row.names=F);
  
  combatout = ComBat(expression_xls='datForComBat.txt', 
                     sample_info_file='infoForComBat.txt', 
                     filter=F, write=F, skip=1,
                     ...);
  return(combatout[,-1]);
}

# sloooooow
.getZerosAcrossConditions = function (dat, conditions) {
  cCounts = table(conditions);
  zDF = data.frame(zCounts=rep(0,nrow(dat)), 
                   sameCond=rep(FALSE,nrow(dat)),
                   sameCondAll=rep(FALSE,nrow(dat)),
                   nameCond=rep(NA,nrow(dat))
                   );
  rownames(zDF) = rownames(dat);
  for (row in 1:nrow(dat)) {
    if (row %% 1000 == 0) {cat(paste0(row,'...'))}
    zInds = which(dat[row, ] == 0);
    if (length(zInds) > 0) {
      zDF$zCounts[row] = length(zInds);
      cDist = table(conditions[zInds]);
      if (length(cDist) == 1) {
        thisCond = names(cDist);
        zDF$sameCond[row] = TRUE;
        zDF$nameCond[row] = thisCond;
        if (cDist == cCounts[match(thisCond, names(cCounts))]) {
          zDF$sameCondAll[row] = TRUE;
        }
      }
    }
  }
  return(zDF);
}