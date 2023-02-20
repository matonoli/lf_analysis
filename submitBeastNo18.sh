#!/bin/bash

nev=1234567890
date=230220
pid=1

#u can force threads by replacing aliroot with 'taskset --cpu-list ${pid} aliroot'

runlist=2016_LHC16k00
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF
cd ..

runlist=2016_LHC16k01
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2016_LHC16k
mkdir workdirMC$runlist
cd workdirMC$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_MC_pass2.list","beastMC_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF



cd ..
runlist=2017_LHC17h00
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..

runlist=2017_LHC17h01
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17h
mkdir workdirMC$runlist
cd workdirMC$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_MC_pass2.list","beastMC_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17i
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17i
mkdir workdirMC$runlist
cd workdirMC$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_MC_pass2.list","beastMC_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17j
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17j
mkdir workdirMC$runlist
cd workdirMC$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_MC_pass2.list","beastMC_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17k00
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..

runlist=2017_LHC17k01
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17k
mkdir workdirMC$runlist
cd workdirMC$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_MC_pass2.list","beastMC_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17l
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17l
mkdir workdirMC$runlist
cd workdirMC$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_MC_pass2.list","beastMC_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17m00
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..

runlist=2017_LHC17m01
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17m
mkdir workdirMC$runlist
cd workdirMC$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_MC_pass2.list","beastMC_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17o00
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..

runlist=2017_LHC17o01
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17o
mkdir workdirMC$runlist
cd workdirMC$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_MC_pass2.list","beastMC_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17r
mkdir workdir$runlist
cd workdir$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_pass2.list","beast_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF

cd ..
runlist=2017_LHC17r
mkdir workdirMC$runlist
cd workdirMC$runlist
. ../prep.sh 09 1
pid=$((pid+1))
aliroot -b -l <<EOF &> ${runlist}.log &
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
root->LoadMacro("MyAnalysisV0unfold.cxx+");
root->LoadMacro("MyAnalysisV0syst.cxx+");
root->LoadMacro("MyAnalysisV0plot.cxx+");
.x ../doAnalysisV0.C(${nev},"0","../lists/beast_${runlist}_MC_pass2.list","beastMC_${date}_${runlist}.root","../rootOutputs/pp16kP2MC_200512_AnalysisResults_hist.root");
.q
EOF