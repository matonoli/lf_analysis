#!/bin/bash

nev=1234567890
date=210604
pid=10

runlist=2018_LHC18b
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
taskset --cpu-list ${pid} aliroot -b -l <<EOF &> ${runlist}.log &
gSystem->Load("libPWGLFspectra");
TInterpreter* root = gInterpreter;
root->LoadMacro("TransverseSpherocity.cxx+");
root->LoadMacro("MyAnalysis.cxx+");
root->LoadMacro("MyHandler.cxx+");
root->LoadMacro("MyEvent.cxx+");
root->LoadMacro("MyTrack.cxx+");
root->LoadMacro("MyParticle.cxx+");
root->LoadMacro("MyV0.cxx++");
root->LoadMacro("MyAnalysisV0.cxx++");
root->LoadMacro("MyAnalysisV0extract.cxx+");
root->LoadMacro("MyAnalysisV0correct.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF
cd ..
runlist=2018_LHC18d
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
taskset --cpu-list ${pid} aliroot -b -l <<EOF &> ${runlist}.log &
gSystem->Load("libPWGLFspectra");
TInterpreter* root = gInterpreter;
root->LoadMacro("TransverseSpherocity.cxx+");
root->LoadMacro("MyAnalysis.cxx+");
root->LoadMacro("MyHandler.cxx+");
root->LoadMacro("MyEvent.cxx+");
root->LoadMacro("MyTrack.cxx+");
root->LoadMacro("MyParticle.cxx+");
root->LoadMacro("MyV0.cxx++");
root->LoadMacro("MyAnalysisV0.cxx++");
root->LoadMacro("MyAnalysisV0extract.cxx+");
root->LoadMacro("MyAnalysisV0correct.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF
cd ..
runlist=2018_LHC18e
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
taskset --cpu-list ${pid} aliroot -b -l <<EOF &> ${runlist}.log &
gSystem->Load("libPWGLFspectra");
TInterpreter* root = gInterpreter;
root->LoadMacro("TransverseSpherocity.cxx+");
root->LoadMacro("MyAnalysis.cxx+");
root->LoadMacro("MyHandler.cxx+");
root->LoadMacro("MyEvent.cxx+");
root->LoadMacro("MyTrack.cxx+");
root->LoadMacro("MyParticle.cxx+");
root->LoadMacro("MyV0.cxx++");
root->LoadMacro("MyAnalysisV0.cxx++");
root->LoadMacro("MyAnalysisV0extract.cxx+");
root->LoadMacro("MyAnalysisV0correct.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF
cd ..
runlist=2018_LHC18f
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
taskset --cpu-list ${pid} aliroot -b -l <<EOF &> ${runlist}.log &
gSystem->Load("libPWGLFspectra");
TInterpreter* root = gInterpreter;
root->LoadMacro("TransverseSpherocity.cxx+");
root->LoadMacro("MyAnalysis.cxx+");
root->LoadMacro("MyHandler.cxx+");
root->LoadMacro("MyEvent.cxx+");
root->LoadMacro("MyTrack.cxx+");
root->LoadMacro("MyParticle.cxx+");
root->LoadMacro("MyV0.cxx++");
root->LoadMacro("MyAnalysisV0.cxx++");
root->LoadMacro("MyAnalysisV0extract.cxx+");
root->LoadMacro("MyAnalysisV0correct.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF
cd ..
runlist=2018_LHC18g
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
taskset --cpu-list ${pid} aliroot -b -l <<EOF &> ${runlist}.log &
gSystem->Load("libPWGLFspectra");
TInterpreter* root = gInterpreter;
root->LoadMacro("TransverseSpherocity.cxx+");
root->LoadMacro("MyAnalysis.cxx+");
root->LoadMacro("MyHandler.cxx+");
root->LoadMacro("MyEvent.cxx+");
root->LoadMacro("MyTrack.cxx+");
root->LoadMacro("MyParticle.cxx+");
root->LoadMacro("MyV0.cxx++");
root->LoadMacro("MyAnalysisV0.cxx++");
root->LoadMacro("MyAnalysisV0extract.cxx+");
root->LoadMacro("MyAnalysisV0correct.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF
cd ..
runlist=2018_LHC18l
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
taskset --cpu-list ${pid} aliroot -b -l <<EOF &> ${runlist}.log &
gSystem->Load("libPWGLFspectra");
TInterpreter* root = gInterpreter;
root->LoadMacro("TransverseSpherocity.cxx+");
root->LoadMacro("MyAnalysis.cxx+");
root->LoadMacro("MyHandler.cxx+");
root->LoadMacro("MyEvent.cxx+");
root->LoadMacro("MyTrack.cxx+");
root->LoadMacro("MyParticle.cxx+");
root->LoadMacro("MyV0.cxx++");
root->LoadMacro("MyAnalysisV0.cxx++");
root->LoadMacro("MyAnalysisV0extract.cxx+");
root->LoadMacro("MyAnalysisV0correct.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF
cd ..
runlist=2018_LHC18m
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
taskset --cpu-list ${pid} aliroot -b -l <<EOF &> ${runlist}.log &
gSystem->Load("libPWGLFspectra");
TInterpreter* root = gInterpreter;
root->LoadMacro("TransverseSpherocity.cxx+");
root->LoadMacro("MyAnalysis.cxx+");
root->LoadMacro("MyHandler.cxx+");
root->LoadMacro("MyEvent.cxx+");
root->LoadMacro("MyTrack.cxx+");
root->LoadMacro("MyParticle.cxx+");
root->LoadMacro("MyV0.cxx++");
root->LoadMacro("MyAnalysisV0.cxx++");
root->LoadMacro("MyAnalysisV0extract.cxx+");
root->LoadMacro("MyAnalysisV0correct.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF
cd ..
runlist=2018_LHC18n
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
taskset --cpu-list ${pid} aliroot -b -l <<EOF &> ${runlist}.log &
gSystem->Load("libPWGLFspectra");
TInterpreter* root = gInterpreter;
root->LoadMacro("TransverseSpherocity.cxx+");
root->LoadMacro("MyAnalysis.cxx+");
root->LoadMacro("MyHandler.cxx+");
root->LoadMacro("MyEvent.cxx+");
root->LoadMacro("MyTrack.cxx+");
root->LoadMacro("MyParticle.cxx+");
root->LoadMacro("MyV0.cxx++");
root->LoadMacro("MyAnalysisV0.cxx++");
root->LoadMacro("MyAnalysisV0extract.cxx+");
root->LoadMacro("MyAnalysisV0correct.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF
cd ..
runlist=2018_LHC18o
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
taskset --cpu-list ${pid} aliroot -b -l <<EOF &> ${runlist}.log &
gSystem->Load("libPWGLFspectra");
TInterpreter* root = gInterpreter;
root->LoadMacro("TransverseSpherocity.cxx+");
root->LoadMacro("MyAnalysis.cxx+");
root->LoadMacro("MyHandler.cxx+");
root->LoadMacro("MyEvent.cxx+");
root->LoadMacro("MyTrack.cxx+");
root->LoadMacro("MyParticle.cxx+");
root->LoadMacro("MyV0.cxx++");
root->LoadMacro("MyAnalysisV0.cxx++");
root->LoadMacro("MyAnalysisV0extract.cxx+");
root->LoadMacro("MyAnalysisV0correct.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF
cd ..
runlist=2018_LHC18p
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
taskset --cpu-list ${pid} aliroot -b -l <<EOF &> ${runlist}.log &
gSystem->Load("libPWGLFspectra");
TInterpreter* root = gInterpreter;
root->LoadMacro("TransverseSpherocity.cxx+");
root->LoadMacro("MyAnalysis.cxx+");
root->LoadMacro("MyHandler.cxx+");
root->LoadMacro("MyEvent.cxx+");
root->LoadMacro("MyTrack.cxx+");
root->LoadMacro("MyParticle.cxx+");
root->LoadMacro("MyV0.cxx++");
root->LoadMacro("MyAnalysisV0.cxx++");
root->LoadMacro("MyAnalysisV0extract.cxx+");
root->LoadMacro("MyAnalysisV0correct.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF
cd ..
runlist=2016_LHC16k
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
taskset --cpu-list ${pid} aliroot -b -l <<EOF &> ${runlist}.log &
gSystem->Load("libPWGLFspectra");
TInterpreter* root = gInterpreter;
root->LoadMacro("TransverseSpherocity.cxx+");
root->LoadMacro("MyAnalysis.cxx+");
root->LoadMacro("MyHandler.cxx+");
root->LoadMacro("MyEvent.cxx+");
root->LoadMacro("MyTrack.cxx+");
root->LoadMacro("MyParticle.cxx+");
root->LoadMacro("MyV0.cxx++");
root->LoadMacro("MyAnalysisV0.cxx++");
root->LoadMacro("MyAnalysisV0extract.cxx+");
root->LoadMacro("MyAnalysisV0correct.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

