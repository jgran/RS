#!/bin/bash

make -j12

#INDIR=/nfs-6/userdata/mt2/V00-00-08_skim_nj2_ht450_mt2gt200_Zinv/
#INDIR=/nfs-6/userdata/mt2/V00-00-11_skim_nj2_ht450_met30_mt2gt200_Zinv/
#INDIR=/nfs-6/userdata/mt2/V00-00-12/
#INDIR=/nfs-6/userdata/mt2/RebalanceAndSmear_V00-00-12_jetpt10/
#INDIR=/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceBabyMaker/
INDIR=/nfs-6/userdata/mt2/RebalanceAndSmear_V00-00-12_jetpt10_genjets_noJEC_V2/

#OUTDIR=output/qcd_undermeasured_mt2_150/
#OUTDIR=output/qcd_undermeasured_mt2_50/
#OUTDIR=output/qcd_undermeasured_mt2_100_ht_500/
#OUTDIR=output/minuit_test/
#OUTDIR=output/100k/
#OUTDIR=output/100k_v2/
#OUTDIR=output/all/
#OUTDIR=output/test/
#OUTDIR=output/jetpt10/
#OUTDIR=output/jetpt10_PUid/
#OUTDIR=output/jetpt10_PUid_alljets/
#OUTDIR=output/jetpt10_PUid_200k/
#OUTDIR=output/additional_plots/
#OUTDIR=output/ht450to1000/
#OUTDIR=output/ht450to1000_rebMET20/
#OUTDIR=output/ht450to1000_rebMET50/
#OUTDIR=output/test_simplex/
#OUTDIR=output/test_gaussian/
#OUTDIR=output/test_fit/
#OUTDIR=output/mht_minimization/
#OUTDIR=output/test_sigmasoft5_50percent/
#OUTDIR=output/test_lowmet/
#OUTDIR=output/test_sigmasoft20_10percent/
#OUTDIR=output/test_sigmasoft5_2percent/
#OUTDIR=output/test_sigmasoft30_2percent/
#OUTDIR=output/test_sigmasoft10_2percent/
#OUTDIR=output/test_sigmasoft15_2percent/
#OUTDIR=output/test_sigmasoft_nvert_2percent/
#OUTDIR=output/test_sigmasoft_nvert_20percent/
#OUTDIR=output/test_sigmasoft_sqrt_nvert_10percent/
#OUTDIR=output/test_sigmasoft20_10percent_ht1000/
#OUTDIR=output/test_sigmasoft20_100percent_ht1000/
#OUTDIR=output/test_sigmasoft20_100percent_ht450to1000/
#OUTDIR=output/test_sigmasoft20_100percent_baseline/
#OUTDIR=output/test_sigmasoft15_100percent_ht450to1000/
OUTDIR=output/test_sigmasoft25_100percent_ht450to1000/

#OUTDIR=output/njets2_ht450to1000_5percent/
#OUTDIR=output/njets2_ht1000_5percent/
#OUTDIR=output/njets3_ht450to1000_5percent/
#OUTDIR=output/njets3_ht1000_5percent/

#OUTDIR=output/rebalance_plots_10percent/

#OUTDIR=output/njets2_ht450to1000_diffmetmht_10percent/
#OUTDIR=output/njets2_ht450to1000_diffmetmht_5percent/

#OUTDIR=output/rebalance_plots_sigmasoft20/

#declare -a Samples=(qcd_pt)
#declare -a Samples=(qcd_pt1000to1400 qcd_pt1400to1800 qcd_pt170to300 qcd_pt1800to2400 qcd_pt2400to3200 qcd_pt300to470 qcd_pt3200 qcd_pt470to600 qcd_pt600to800 qcd_pt800to1000 qcd_pt120to170)
#declare -a Samples=(qcd_pt1800to2400 qcd_pt300to470)
#declare -a Samples=(qcd_pt600to800)
#declare -a Samples=(qcd_pt3200)
#declare -a Samples=(qcd_pt170to300)
#declare -a Samples=(qcd_pt1800to2400)
#declare -a Samples=(sntMT2Baby)
#declare -a Samples=(qcd_pt120to170)
#declare -a Samples=(qcd_pt300to470 qcd_pt3200 qcd_pt470to600 qcd_pt600to800 qcd_pt800to1000)
#declare -a Samples=(qcd_pt1000to1400 qcd_pt1400to1800 qcd_pt170to300 qcd_pt1800to2400 qcd_pt2400to3200)

#declare -a Samples=(qcd_pt120to170 qcd_pt170to300 qcd_pt300to470 qcd_pt470to600 qcd_pt600to800 qcd_pt800to1000 )
#declare -a Samples=(qcd_pt1000to1400 qcd_pt1400to1800 qcd_pt1800to2400 qcd_pt2400to3200 qcd_pt3200 )
declare -a Samples=(qcd_pt120to170 qcd_pt170to300)

mkdir -p ${OUTDIR}

for SAMPLE in ${Samples[@]};
  do echo root -b -q -l doAll.C\(\"${INDIR}\",\"${SAMPLE}\",\"${OUTDIR}\"\)
  nohup nice -5 root -b -q -l doAll.C\(\"${INDIR}\",\"${SAMPLE}\",\"${OUTDIR}\"\) >& log_${SAMPLE}.txt &
done
