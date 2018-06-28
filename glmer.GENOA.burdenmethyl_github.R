####################################
#	Using Genetic Burden Scores for Gene-by-Methylation Interaction Analysis
#	on Metabolic Syndrome in African Americans
#
#	+Jacquelyn Y. Taylor, Rory Meyers College of Nursing, jt139@nyu.edu
#	+Erin B. Ware, PhD, MPH, ebakshis@umich.edu
#	Michelle L. Wright, PhD, RN, mlwrig5@emory.edu
#	Jennifer A. Smith, PhD, MPH, smjenn@umich.edu
# 	Sharon L.R. Kardia, PhD, skardia@umich.edu
#
#	+Co-first authors
####################################


library(lme4)

########BURDEN SCORES###########
########BURDEN SCORES###########
########BURDEN SCORES###########


#READ IN PHENOTYPE FILE
#(example:     data.pheno<-read.table('/directory/phenotype.csv', sep=",", header=T,stringsAsFactors=F)     )
data.pheno<-read.csv("pheno.geno.methyl.739.csv")
data.geno<-read.csv("AA.burdenscores.csv")
lme_data<-merge(data.pheno, data.geno, by="PIN")

ncol.pheno <- 127

#SET OUTPUT DIRECTORY
home_dir <- './'

ncol.lme_data <- ncol(lme_data)

#DEFINE SNP COLUMNS
snps <- colnames(lme_data)[110:127]

#LIST OUT N-1 BURDEN SCORE VALUES, HOLD ONE OUT TO SET FILE
#"CA10_05cont"
snps<-c("KSR2_05cont", "LPL_05cont","APOA5_05cont","LCAT_05cont", "ABCG8_05cont","DHODH_05cont","CTNNA3_01cont","RALYL_05cont")

#DEFINE COVARIATES, SET OUTCOME AND CLUSTER VARIABLE
group_vars=c('PC1','PC2','PC3','PC4','exf_age', 'GENDER_n')
random="pedid"
measure_var="mets"

#CREATE CASE 1
snp_id="CA10_05cont" 
	cols <- c(random, group_vars, measure_var, snp_id)
	lme_data2 <- na.omit(lme_data[, cols])
	if(is.null(group_vars)) {
		rhs <- snp_id
	} else {
		rhs <- paste(snp_id, '+', paste(group_vars, collapse='+'))
	}
	lme_formula <- as.formula(paste(measure_var, '~', rhs,"+(1|pedid)"))
	lmeResult<-glmer(lme_formula, lme_data2, family=binomial, control=glmerControl(optimizer="bobyqa"))

#CHECK RESULTS FOR CASE 1	
summary(lmeResult)
summary(lmeResult)$coefficients

#SET UP RESULTS MATRIX
result<-matrix(NA,1,10)
colnames(result)<-c('burden','n','int_beta', 'int_beta_se', 'int_beta_Z', 'int_beta_p','burden_beta', 'burden_beta_se', 'burden_beta_Z', 'burden_beta_p')
result[1,1]=snp_id
result[1,2]=dim(lme_data2)[1]
result[1,3]=summary(lmeResult)$coefficients[1,1]
result[1,4]=summary(lmeResult)$coefficients[1,2]
result[1,5]=summary(lmeResult)$coefficients[1,3]
result[1,6]=summary(lmeResult)$coefficients[1,4]
result[1,7]=summary(lmeResult)$coefficients[2,1]
result[1,8]=summary(lmeResult)$coefficients[2,2]
result[1,9]=summary(lmeResult)$coefficients[2,3]
result[1,10]=summary(lmeResult)$coefficients[2,4]
result

#LOOP THROUGH ALL SNPS
for (snp_id in snps){
	cols <- c(random, group_vars, measure_var, snp_id)
	lme_data2 <- na.omit(lme_data[, cols])

	rhs <- paste(snp_id, '+', paste(group_vars, collapse='+'))

	lme_formula <- as.formula(paste(measure_var, '~', rhs,"+(1|pedid)"))
	
	lmeResult<-glmer(lme_formula, lme_data2, family=binomial, control=glmerControl(optimizer="bobyqa"))
	results<-matrix(NA,1,10)
	results[1,1]=snp_id
	results[1,2]=dim(lme_data2)[1]
	results[1,3]=summary(lmeResult)$coefficients[1,1]
	results[1,4]=summary(lmeResult)$coefficients[1,2]
	results[1,5]=summary(lmeResult)$coefficients[1,3]
	results[1,6]=summary(lmeResult)$coefficients[1,4]
	results[1,7]=summary(lmeResult)$coefficients[2,1]
	results[1,8]=summary(lmeResult)$coefficients[2,2]
	results[1,9]=summary(lmeResult)$coefficients[2,3]
	results[1,10]=summary(lmeResult)$coefficients[2,4]

	result<-rbind(result, results)
	print(paste0("done with gene: ",snp_id))
}

dim(result)
result

#SET UP OUTPUT FILE
date <- format(Sys.time(), "%Y%m%d")
out_dir <- paste(home_dir, 'Results', sep='/')
out_fn_pre <- paste(out_dir, 'genoa.contburden.mets', sep='/') 
out_fn_suf <- paste(date, 'csv', sep='.')
cov_str <- paste(group_vars, collapse='_')
output_fn = paste(out_fn_pre, cov_str, out_fn_suf, sep='.')

#Write output table
write.table(result, output_fn, row.names=FALSE, sep=',', quote=FALSE)


########METHYLATION###########
########METHYLATION###########
########METHYLATION###########

#READ IN YOUR PHENOTYPE FILE
data.pheno<-read.csv("./data/pheno.geno.methyl.739.csv")

#SET OUTPUT DIRECTORY
home_dir <- './'
lme_data<-data.pheno

#LIST OUT N-1 BURDEN SCORE VALUES, HOLD ONE OUT TO SET FILE
#"cg00132141" 
methly<-c("cg09538287","cg08918749","cg22108175","cg02157083","cg25682080","cg01489608","cg26924825","cg07817698","cg22381196")

##DEFINE COVARIATES, OUTCOME, CLUSTER VARIABLE
group_vars=c('exf_age', 'GENDER_n')
random="pedid"
measure_var="mets"

#CREATE CASE 1

M_id="cg00132141" 
	cols <- c(random, group_vars, measure_var, M_id)
	lme_data2 <- na.omit(lme_data[, cols])
	if(is.null(group_vars)) {
		rhs <- M_id
	} else {
		rhs <- paste(M_id, '+', paste(group_vars, collapse='+'))
	}
	lme_formula <- as.formula(paste(measure_var, '~', rhs,"+(1|pedid)"))
	lmeResult<-glmer(lme_formula, lme_data2, family=binomial, control=glmerControl(optimizer="bobyqa"))
#CHECK RESULTS FOR CASE 1
summary(lmeResult)
summary(lmeResult)$coefficients

#SET UP RESULTS MATRIX
result<-matrix(NA,1,10)
colnames(result)<-c('methyl','n','int_beta', 'int_beta_se', 'int_beta_Z', 'int_beta_p','methyl_beta', 'methyl_beta_se', 'methyl_beta_Z', 'methyl_beta_p')
	result[1,1]=M_id
	result[1,2]=dim(lme_data2)[1]
	result[1,3]=summary(lmeResult)$coefficients[1,1]
	result[1,4]=summary(lmeResult)$coefficients[1,2]
	result[1,5]=summary(lmeResult)$coefficients[1,3]
	result[1,6]=summary(lmeResult)$coefficients[1,4]
	result[1,7]=summary(lmeResult)$coefficients[2,1]
	result[1,8]=summary(lmeResult)$coefficients[2,2]
	result[1,9]=summary(lmeResult)$coefficients[2,3]
	result[1,10]=summary(lmeResult)$coefficients[2,4]
	result

#LOOP THROUGH ALL M* SITES
for (M_id in methly){
	cols <- c(random, group_vars, measure_var, M_id)
	lme_data2 <- na.omit(lme_data[, cols])

	rhs <- paste(M_id, '+', paste(group_vars, collapse='+'))

	lme_formula <- as.formula(paste(measure_var, '~', rhs,"+(1|pedid)"))
	
	lmeResult<-glmer(lme_formula, lme_data2, family=binomial, control=glmerControl(optimizer="bobyqa"))
	results<-matrix(NA,1,10)
	results[1,1]=M_id
	results[1,2]=dim(lme_data2)[1]
	results[1,3]=summary(lmeResult)$coefficients[1,1]
	results[1,4]=summary(lmeResult)$coefficients[1,2]
	results[1,5]=summary(lmeResult)$coefficients[1,3]
	results[1,6]=summary(lmeResult)$coefficients[1,4]
	results[1,7]=summary(lmeResult)$coefficients[2,1]
	results[1,8]=summary(lmeResult)$coefficients[2,2]
	results[1,9]=summary(lmeResult)$coefficients[2,3]
	results[1,10]=summary(lmeResult)$coefficients[2,4]

	result<-rbind(result, results)
	print(paste0("done with site: ",M_id))
}

dim(result)
result

#SET UP OUTPUT FILE
date <- format(Sys.time(), "%Y%m%d")
out_dir <- paste(home_dir, 'Results', sep='/')
out_fn_pre <- paste(out_dir, 'genoa.methyl.mets', sep='/') #consider replacing the word phenotype with YOUR phenotype
out_fn_suf <- paste(date, 'csv', sep='.')

cov_str <- paste(group_vars, collapse='_')
output_fn = paste(out_fn_pre, cov_str, out_fn_suf, sep='.')

#WRITE OUTPUT FILE
write.table(result, output_fn, row.names=FALSE, sep=',', quote=FALSE)




########INTERACTION###########
########INTERACTION###########
########INTERACTION###########

#		#M* SITES WITH P<0.2
#		CHR CpG			GENE
#		16	cg22381196	DHODH
#		10	cg00132141	CTNNA3
#		16	cg07817698	DHODH
#		16	cg01489608	LCAT
#		11	cg25682080	APOA5
#		
#		#LCAT_05cont is the only burden score with P<0.2




summary(lme_data[,c("LCAT_05cont","cg01489608")])


lmeResult1a <- glmer(mets~PC1+PC2+PC3+PC4+exf_age+GENDER_n+LCAT_05cont+cg01489608+(1|pedid),lme_data, family=binomial, control=glmerControl(optimizer="bobyqa")) 
lmeResult1b <- glmer(mets~PC1+PC2+PC3+PC4+exf_age+GENDER_n+LCAT_05cont+cg01489608+LCAT_05cont*cg01489608+(1|pedid),lme_data, family=binomial, control=glmerControl(optimizer="bobyqa")) 

results<-matrix(NA,2,19)

colnames(results)<-c('methyl','gene','n','int_beta', 'int_beta_se', 'int_beta_Z', 'int_beta_p','gene_beta', 'gene_beta_se', 'gene_beta_Z', 'gene_beta_p','methyl_beta', 'methyl_beta_se', 'methyl_beta_Z', 'methyl_beta_p','GxE_beta', 'GxE_beta_se', 'GxE_beta_Z', 'GxE_beta_p')

	results[1,1]="cg01489608"
	results[1,2]="LCAT_05cont"
	results[1,3]=dim(lme_data2)[1]
	results[1,4]=summary(lmeResult1a)$coefficients[1,1]
	results[1,5]=summary(lmeResult1a)$coefficients[1,2]
	results[1,6]=summary(lmeResult1a)$coefficients[1,3]
	results[1,7]=summary(lmeResult1a)$coefficients[1,4]
	results[1,8]= summary(lmeResult1a)$coefficients[8,1]
	results[1,9]= summary(lmeResult1a)$coefficients[8,2]
	results[1,10]=summary(lmeResult1a)$coefficients[8,3]
	results[1,11]=summary(lmeResult1a)$coefficients[8,4]
	results[1,12]=summary(lmeResult1a)$coefficients[9,1]
	results[1,13]=summary(lmeResult1a)$coefficients[9,2]
	results[1,14]=summary(lmeResult1a)$coefficients[9,3]
	results[1,15]=summary(lmeResult1a)$coefficients[9,4]
	results[1,16]=NA
	results[1,17]=NA
	results[1,18]=NA
	results[1,19]=NA

	results[2,1]="cg01489608"
	results[2,2]="LCAT_05cont"
	results[2,3]=dim(lme_data2)[1]
	results[2,4]=summary(lmeResult1b)$coefficients[1,1]
	results[2,5]=summary(lmeResult1b)$coefficients[1,2]
	results[2,6]=summary(lmeResult1b)$coefficients[1,3]
	results[2,7]=summary(lmeResult1b)$coefficients[1,4]
	results[2,8]= summary(lmeResult1b)$coefficients[8,1]
	results[2,9]= summary(lmeResult1b)$coefficients[8,2]
	results[2,10]=summary(lmeResult1b)$coefficients[8,3]
	results[2,11]=summary(lmeResult1b)$coefficients[8,4]
	results[2,12]=summary(lmeResult1b)$coefficients[9,1]
	results[2,13]=summary(lmeResult1b)$coefficients[9,2]
	results[2,14]=summary(lmeResult1b)$coefficients[9,3]
	results[2,15]=summary(lmeResult1b)$coefficients[9,4]
	results[2,16]=summary(lmeResult1b)$coefficients[10,1]
	results[2,17]=summary(lmeResult1b)$coefficients[10,2]
	results[2,18]=summary(lmeResult1b)$coefficients[10,3]
	results[2,19]=summary(lmeResult1b)$coefficients[10,4]

results
	
#SET UP OUTPUT FILE
date <- format(Sys.time(), "%Y%m%d")
out_dir <- paste(home_dir, 'Results', sep='/')
out_fn_pre <- paste(out_dir, 'genoa.continteraction.mets', sep='/') #consider replacing the word phenotype with YOUR phenotype
out_fn_suf <- paste(date, 'csv', sep='.')

cov_str <- paste(group_vars, collapse='_')
output_fn = paste(out_fn_pre, cov_str, out_fn_suf, sep='.')

#WRITE OUTPUT FILE
write.table(results, output_fn, row.names=FALSE, sep=',', quote=FALSE)
