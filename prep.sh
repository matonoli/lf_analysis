#!/bin/sh
# PARAMETERS
CUT="cuts"$1".h"
MODE=$2
TREES=1
ESDS=2

# COPY FILES
cp ../TransverseSpherocity/TransverseSpherocity.cxx .
cp ../TransverseSpherocity/TransverseSpherocity.h .		# extra classes

cp ../MyKit/MyAnalysis.cxx .
cp ../MyKit/MyAnalysis.h .								# MyKit
cp ../MyKit/MyHandler.cxx .
cp ../MyKit/MyHandler.h .
cp ../MyKit/MyEvent.cxx .
cp ../MyKit/MyEvent.h .
cp ../MyKit/MyTrack.cxx .
cp ../MyKit/MyTrack.h .
cp ../MyKit/MyParticle.cxx .
cp ../MyKit/MyParticle.h .
cp ../MyKit/MyV0.cxx .
cp ../MyKit/MyV0.h .
cp ../MyKit/Analyses/MyAnalysisV0.cxx .
cp ../MyKit/Analyses/MyAnalysisV0.h .					# MyKit analyses
cp ../MyKit/Analyses/MyAnalysisV0extract.cxx .
cp ../MyKit/Analyses/MyAnalysisV0extract.h .
cp ../MyKit/Analyses/MyAnalysisV0correct.cxx .
cp ../MyKit/Analyses/MyAnalysisV0correct.h .
cp ../MyKit/Analyses/MyAnalysisV0syst.cxx .
cp ../MyKit/Analyses/MyAnalysisV0syst.h .
cp ../MyKit/Analyses/MyAnalysisV0plot.cxx .
cp ../MyKit/Analyses/MyAnalysisV0plot.h .
cp ../AliAnalysisTaskMyTask.cxx .
cp ../AliAnalysisTaskMyTask.h .							# alice task
cp ../AddMyTask.C .

cp ../$CUT .											# cut file

# DEFINE INSTRUCTIONS FOR COMPILING
echo "#ifndef COMPINSTRUCTIONS_H" > compInstructions.h
echo "#define COMPINSTRUCTIONS_H" >> compInstructions.h
echo "#define INPUTFORMAT $MODE" >> compInstructions.h
echo "// $TREES : local aurora trees" >> compInstructions.h
echo "// $ESDS : esd files \n" >> compInstructions.h
echo "#include \"$CUT\"" >> compInstructions.h
echo '#endif' >> compInstructions.h