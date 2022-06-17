### snippets for LDSC and GenomicSEM
# ldsc snippet1 
module load plink
plink --bfile ./data/round8_unpruned --cm-map ./ldsc/HS_genetic_map.SexAvg.chr@.SHAPEIT.txt --make-bed --out ./ldsc/round8_unprunedGMap
# ldsc snippet2
module load plink
for i in {1..20}
do
   plink --bfile ./ldsc/round8_unprunedGMap --chr $i --make-bed --out ./ldsc/gMapChr/round8_unprunedGMapChr$i
done
# ldsc snippet3
source activate ldsc
python ./software/ldsc/ldsc.py --bfile ./ldsc/round8_unprunedGMap --l2 --ld-wind-cm 1 --out ./ldsc/ldscPy/round8_unprunedLdscPyWd1cm
python ./software/ldsc/ldsc.py --bfile ./ldsc/round8_unprunedGMap --l2 --ld-wind-cm 5 --out ./ldsc/ldscPy/round8_unprunedLdscPyWd5cm
python ./software/ldsc/ldsc.py --bfile ./ldsc/round8_unprunedGMap --l2 --ld-wind-cm 10 --out ./ldsc/ldscPy/round8_unprunedLdscPyWd10cm
python ./software/ldsc/ldsc.py --bfile ./ldsc/round8_unprunedGMap --l2 --ld-wind-cm 20 --out ./ldsc/ldscPy/round8_unprunedLdscPyWd20cm
# ldsc snippet4
source activate ldsc
python ./software/ldsc/ldsc.py --bfile ./ldsc/round8_unprunedGMap --l2 --ld-wind-cm 1 --out ./ldsc/ldscPy/round8_unprunedLdscPy
source deactivate
# ldsc snippet5
source activate ldsc
for i in {1..20}
do
   python ./software/ldsc/ldsc.py --bfile ./ldsc/gMapChr/round8_unprunedGMapChr$i --l2 --ld-wind-cm 1 --out ./ldsc/ldscChr/round8_unprunedGMapLdscChr$i
done
gzip -k ldscChr
# ldsc snippet6
source activate ldsc
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d1_total_distResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d1_total_dist
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d1_avg_velResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d1_avg_vel
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d2_total_distResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d2_total_dist
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d2_avg_velResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d2_avg_vel
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d3_total_distResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d3_total_dist
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d3_avg_velResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d3_avg_vel
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d3_total_hw_boutsResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d3_total_hw_bouts
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d3_total_hw_durResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d3_total_hw_dur
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d7_total_distResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d7_total_dist
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d7_avg_velResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d7_avg_vel
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d7_total_hw_boutsResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d7_total_hw_bouts
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d7_total_hw_durResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d7_total_hw_dur
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d8_total_distResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d8_total_dist
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_d8_avg_velResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_d8_avg_vel
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_degree_sens_distResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_degree_sens_dist
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_degree_sens_stereotypyResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_degree_sens_stereotypy
python ./software/ldsc/munge_sumstats.py --sumstats ./data/cccResiduals_con_condi_dist_diff_scoreResid.loco.mlma --N 1376 --out  ./ldsc/sumstats/cccResiduals_con_condi_dist_diff_score
# ldsc snippet7
source activate ldsc
python ./software/ldsc/ldsc.py --h2 ./ldsc/sumstats/cccResiduals_con_condi_dist_diff_score.sumstats.gz --ref-ld ./ldsc/ldscPy/round8_unprunedLdscPy --w-ld ./ldsc/ldscPy/round8_unprunedLdscPy --out ./ldsc/ldscPy/cccResiduals_con_condi_dist_diff_scoreLdscpyh2
# GenomicSEM snippet1
traits <- c("./ldsc/sumstats/cccResiduals_d1_avg_vel.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_d2_avg_vel.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_d3_avg_vel.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_d7_avg_vel.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_d8_avg_vel.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_d1_total_dist.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_d2_total_dist.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_d3_total_dist.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_d7_total_dist.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_d8_total_dist.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_con_condi_dist_diff_score.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_degree_sens_dist.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_d3_total_hw_bouts.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_d3_total_hw_dur.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_d7_total_hw_bouts.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_d7_total_hw_dur.sumstats.gz",
            "./ldsc/sumstats/cccResiduals_degree_sens_stereotypy.sumstats.gz")
sample.prev <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
population.prev <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
ld <-"./ldsc/ldscChr/round8_unprunedGMapLdscChr"
wld <-"./ldsc/ldscChr/round8_unprunedGMapLdscChr"
trait.names<-c("d1Vel", "d2Vel", "d3Vel", "d7Vel", "d8Vel", "d1Dist", "d2Dist", "d3Dist", "d7Dist",
               "d8Dist", "CCDDS", "DSD", "d3Bouts", "d3Dur", "d7Bouts", "d7Dur", "DSS") 
cccLDSCoutputGSEM1 <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, chr=20) 
save(cccLDSCoutputGSEM1, file="./ldsc/cccLDSCoutputGSEM1.RData") ## check
# GenomicSEM snippet2
load("./ldsc/cccLDSCoutputGSEM1.RData")
modelToTest <-'
Vel =~ d1Vel + d2Vel + d3Vel + d7Vel + d8Vel
Dist =~ d1Dist + d2Dist + d3Dist + d7Dist + d8Dist + CCDDS + DSD
Act =~ d3Bouts + d3Dur + d7Bouts + d7Dur + DSS
Vel ~~ Dist
Vel ~~ Act
Dist ~~ Act
'
cccResultsGSEM2<-usermodel(cccLDSCoutputGSEM1, estimation = "DWLS", model = modelToTest, CFIcalc = TRUE,
                           std.lv = TRUE, imp_cov = FALSE)#smoothCheck 
save(cccResultsGSEM2, file="./ldsc/cccResultsGSEM2.RData") 
# GenomicSEM snippet3
munge(c("./data/cccResiduals_d1_avg_velResid.loco.mlma",
        "./data/cccResiduals_d2_avg_velResid.loco.mlma",
        "./data/cccResiduals_d3_avg_velResid.loco.mlma",
        "./data/cccResiduals_d7_avg_velResid.loco.mlma",
        "./data/cccResiduals_d8_avg_velResid.loco.mlma",
        "./data/cccResiduals_d1_total_distResid.loco.mlma",
        "./data/cccResiduals_d2_total_distResid.loco.mlma",
        "./data/cccResiduals_d3_total_distResid.loco.mlma",
        "./data/cccResiduals_d8_total_distResid.loco.mlma",
        "./data/cccResiduals_con_condi_dist_diff_scoreResid.loco.mlma",
        "./data/cccResiduals_degree_sens_distResid.loco.mlma",
        "./data/cccResiduals_d7_total_hw_boutsResid.loco.mlma",
        "./data/cccResiduals_d7_total_hw_durResid.loco.mlma",
        "./data/cccResiduals_degree_sens_stereotypyResid.loco.mlma"),
      hm3 = "./ldsc/round8_unprunedReferenceFile.txt", 
      trait.names=c("d1Vel", "d2Vel", "d3Vel", "d7Vel", "d8Vel", "d1Dist", "d2Dist", "d3Dist", "d8Dist",
                    "CCDDS", "DSD", "d7Bouts", "d7Dur", "DSS"), 
      c(1376, 1376, 1376, 1376, 1376, 1376, 1376, 1376, 1376, 1376, 1376, 1376, 1376, 1376))
# GenomicSEM snippet4
traits <- c("./ldsc/mungeSEM4/d1Vel.sumstats.gz",
            "./ldsc/mungeSEM4/d2Vel.sumstats.gz",
            "./ldsc/mungeSEM4/d3Vel.sumstats.gz",
            "./ldsc/mungeSEM4/d7Vel.sumstats.gz",
            "./ldsc/mungeSEM4/d8Vel.sumstats.gz",
            "./ldsc/mungeSEM4/d1Dist.sumstats.gz",
            "./ldsc/mungeSEM4/d2Dist.sumstats.gz",
            "./ldsc/mungeSEM4/d3Dist.sumstats.gz",
            "./ldsc/mungeSEM4/d8Dist.sumstats.gz",
            "./ldsc/mungeSEM4/CCDDS.sumstats.gz",
            "./ldsc/mungeSEM4/DSD.sumstats.gz",
            "./ldsc/mungeSEM4/d7Bouts.sumstats.gz",
            "./ldsc/mungeSEM4/d7Dur.sumstats.gz",
            "./ldsc/mungeSEM4/DSS.sumstats.gz")
sample.prev <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
population.prev <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
ld <-"./ldsc/ldscChr/round8_unprunedGMapLdscChr"
wld <-"./ldsc/ldscChr/round8_unprunedGMapLdscChr"
trait.names<-c("d1Vel", "d2Vel", "d3Vel", "d7Vel", "d8Vel", "d1Dist", "d2Dist", "d3Dist", "d8Dist",
               "CCDDS", "DSD", "d7Bouts", "d7Dur", "DSS") 
cccLDSCoutputGSEM5 <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, chr=20) 
save(cccLDSCoutputGSEM5, file="./ldsc/cccLDSCoutputGSEM5.RData") ## check
# GenomicSEM snippet5
se.logit <-c(F,F,F,F,F,F,F,F,F,F,F,F,F,F) 
### sumstats function
files <- list(".data/d1_avg_velResid.loco.mlma",
              ".data/cccResiduals_d2_avg_velResid.loco.mlma",
              ".data/cccResiduals_d3_avg_velResid.loco.mlma",
              ".data/cccResiduals_d7_avg_velResid.loco.mlma",
              ".data/cccResiduals_d8_avg_velResid.loco.mlma",
              ".data/cccResiduals_d1_total_distResid.loco.mlma",
              ".data/cccResiduals_d2_total_distResid.loco.mlma",
              ".data/cccResiduals_d3_total_distResid.loco.mlma",
              ".data/cccResiduals_d8_total_distResid.loco.mlma",
              ".data/cccResiduals_con_condi_dist_diff_scoreResid.loco.mlma",
              ".data/cccResiduals_degree_sens_distResid.loco.mlma",
              ".data/cccResiduals_d7_total_hw_boutsResid.loco.mlma",
              ".data/cccResiduals_d7_total_hw_durResid.loco.mlma",
              ".data/cccResiduals_degree_sens_stereotypyResid.loco.mlma") 
trait.names <- c("d1Vel", "d2Vel", "d3Vel", "d7Vel", "d8Vel", "d1Dist", "d2Dist", "d3Dist", "d8Dist", "CCDDS", "DSD", "d7Bouts", "d7Dur", "DSS")
ref <- "./reference/round8_unprunedReferenceFile.txt"
cccSumstatsSEM6 <- sumstats(files=files,ref=ref, trait.names=trait.names,se.logit=se.logit, maf.filter=0.01) 
save(cccSumstatsSEM6, file="./ldsc/cccSumstatsGSEM6.RData") ## check
# GenomicSEM snippet6
load("./ldsc/cccLDSCoutputGSEM5.RData")
load("./ldsc/cccSumstatsGSEM6.RData")
#Model to test
model <-'
Vel =~ d1Vel + d2Vel + d3Vel + d7Vel + d8Vel
Dist =~ d1Dist + d2Dist + d3Dist + d8Dist + CCDDS + DSD
Act =~ d7Bouts + d7Dur + DSS
Vel ~~ Dist
Dist ~~ Act
Vel ~~ Act
Vel ~ SNP 
Dist ~ SNP
Act ~ SNP'
#run the multivariate GWAS
cccCorrelatedFactorsGSEM7<-userGWAS(covstruc = cccLDSCoutputGSEM5, SNPs = cccSumstatsSEM6, estimation = "DWLS",
                                    model = model, printwarn = TRUE, sub=c("Vel~SNP", "Dist~SNP", "Act~SNP"),
                                    cores = 8, toler = FALSE, SNPSE = FALSE, parallel = FALSE,GC="standard",MPI=FALSE) 
save(cccCorrelatedFactorsGSEM7, file="./ldsc/cccCorrelatedFactorsGSEM7.RData") 
