#!/bin/bash

make -j12

#INDIR=/nfs-6/userdata/mt2/V00-00-08_skim_nj2_ht450_mt2gt200_Zinv/
#INDIR=/nfs-6/userdata/mt2/V00-00-11_skim_nj2_ht450_met30_mt2gt200_Zinv/
#INDIR=/nfs-6/userdata/mt2/RebalanceAndSmear_V00-00-12_jetpt10/
#INDIR=/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceBabyMaker/
INDIR=/nfs-6/userdata/mt2/RebalanceAndSmear_V00-00-12_jetpt10_genjets_noJEC_V2/

#OUTDIR=output/qcd_undermeasured_mt2_150/
#OUTDIR=output/qcd_undermeasured_mt2_50/
#OUTDIR=output/qcd_undermeasured_mt2_100_ht_500/
#OUTDIR=output/minuit_test/
#OUTDIR=output/100k/
#OUTDIR=output/jetpt10/
#OUTDIR=output/jetpt10_PUid/
#OUTDIR=output/jetpt10_PUid_200k/
#OUTDIR=output/jetpt10_PUid_all/
#OUTDIR=output/test/
#OUTDIR=output/test_sigmasoft5_10percent/
#OUTDIR=output/test_sigmasoft5_50percent/
#OUTDIR=output/test_simplex/
#OUTDIR=output/test_gaussian/
#OUTDIR=output/test_fit/
#OUTDIR=output/mht_minimization/
#OUTDIR=output/met_minimization/
#OUTDIR=output/met_minimization_sigmasoft20_5percent/
#OUTDIR=output/met_minimization_sigmasoft20_2percent/
#OUTDIR=output/met_minimization_sigmasoft5_2percent/
#OUTDIR=output/met_minimization_sigmasoft30_2percent/
#OUTDIR=output/met_minimization_sigmasoft10_2percent/
#OUTDIR=output/met_minimization_sigmasoft15_2percent/
#OUTDIR=output/met_minimization_sigmasoft_nvert_2percent/
#OUTDIR=output/met_minimization_sigmasoft_sqrtnvert_10percent/
#OUTDIR=output/met_minimization_sigmasoft20_100percent/
#OUTDIR=output/met_minimization_sigmasoft15_100percent/
OUTDIR=output/met_minimization_sigmasoft25_100percent/

#declare -a Samples=(qcd_pt1000to1400 qcd_pt1400to1800 qcd_pt170to300 qcd_pt1800to2400 qcd_pt2400to3200 qcd_pt300to470 qcd_pt3200 qcd_pt470to600 qcd_pt600to800 qcd_pt800to1000)
#declare -a Samples=(qcd_pt)
#declare -a Samples=(qcd_pt1000to1400)
#declare -a Samples=(qcd_pt600to800)
#declare -a Samples=(qcd_pt3200)
#declare -a Samples=(qcd_pt170to300)
#declare -a Samples=(qcd_pt300to470)
#declare -a Samples=(qcd_pt120to170)
declare -a Samples=(qcd_pt300to470 qcd_pt3200 qcd_pt470to600 qcd_pt600to800 qcd_pt800to1000)
#declare -a Samples=(qcd_pt1000to1400 qcd_pt1400to1800 qcd_pt170to300 qcd_pt1800to2400 qcd_pt2400to3200)

mkdir -p ${OUTDIR}

for SAMPLE in ${Samples[@]};
  do echo root -b -q -l doAll.C\(\"${INDIR}\",\"${SAMPLE}\",\"${OUTDIR}\"\)
  #nohup nice -5 root -b -q -l doAll.C\(\"${INDIR}\",\"${SAMPLE}\",\"${OUTDIR}\"\) >& log_${SAMPLE}.txt &
  nohup nice -5 root -b -q -l doAll.C\(\"${INDIR}\",\"${SAMPLE}\",\"${OUTDIR}\"\) > /dev/null 2>&1 &
  #root -l doAll.C\(\"${INDIR}\",\"${SAMPLE}\",\"${OUTDIR}\"\)
done
