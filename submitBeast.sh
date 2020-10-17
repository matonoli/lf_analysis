#!/bin/bash
echo "What runlist? Like: 2016_LHC16k or 2016_LHC16k_MC"
read runlist
echo "How many events per process?"
read nev
echo "Which process id?"
read pid
echo "What day is it today 201017?"
read date

mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 07 1
taskset --cpu-list ${pid} aliroot -b -l <<EOF > ${runlist}.log
.x ../compile.C;
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2_SHORT.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root")
.q
EOF

