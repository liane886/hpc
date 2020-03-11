##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=111
ConfigurationName      :=Debug
WorkspacePath          :=/home/li/Desktop/hpc
ProjectPath            :=/home/li/Desktop/hpc/111
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=li
Date                   :=11/03/20
CodeLitePath           :=/home/li/.codelite
LinkerName             :=g++
SharedObjectLinkerName :=g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.o.i
DebugSwitch            :=-gstab
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E 
ObjectsFileList        :="111.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := $(LibrarySwitch)blas,lapack,mpi,cblas 
ArLibs                 :=  "blas,lapack,mpi,cblas" 
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := ar rcus
CXX      := g++
CC       := gcc
CXXFLAGS :=  -g -O0 -Wall $(Preprocessors)
CFLAGS   :=  -g -O0 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCmngMPI.cpp$(ObjectSuffix) $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCpoissonSolver.cpp$(ObjectSuffix) $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCprogram_options.cpp$(ObjectSuffix) $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavity.cpp$(ObjectSuffix) $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavitySolver.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


$(IntermediateDirectory)/.d:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCmngMPI.cpp$(ObjectSuffix): ../../LidDrivenCavitySolver-master/src/LDCmngMPI.cpp $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCmngMPI.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/li/Desktop/LidDrivenCavitySolver-master/src/LDCmngMPI.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCmngMPI.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCmngMPI.cpp$(DependSuffix): ../../LidDrivenCavitySolver-master/src/LDCmngMPI.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCmngMPI.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCmngMPI.cpp$(DependSuffix) -MM ../../LidDrivenCavitySolver-master/src/LDCmngMPI.cpp

$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCmngMPI.cpp$(PreprocessSuffix): ../../LidDrivenCavitySolver-master/src/LDCmngMPI.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCmngMPI.cpp$(PreprocessSuffix) ../../LidDrivenCavitySolver-master/src/LDCmngMPI.cpp

$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCpoissonSolver.cpp$(ObjectSuffix): ../../LidDrivenCavitySolver-master/src/LDCpoissonSolver.cpp $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCpoissonSolver.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/li/Desktop/LidDrivenCavitySolver-master/src/LDCpoissonSolver.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCpoissonSolver.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCpoissonSolver.cpp$(DependSuffix): ../../LidDrivenCavitySolver-master/src/LDCpoissonSolver.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCpoissonSolver.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCpoissonSolver.cpp$(DependSuffix) -MM ../../LidDrivenCavitySolver-master/src/LDCpoissonSolver.cpp

$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCpoissonSolver.cpp$(PreprocessSuffix): ../../LidDrivenCavitySolver-master/src/LDCpoissonSolver.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCpoissonSolver.cpp$(PreprocessSuffix) ../../LidDrivenCavitySolver-master/src/LDCpoissonSolver.cpp

$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCprogram_options.cpp$(ObjectSuffix): ../../LidDrivenCavitySolver-master/src/LDCprogram_options.cpp $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCprogram_options.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/li/Desktop/LidDrivenCavitySolver-master/src/LDCprogram_options.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCprogram_options.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCprogram_options.cpp$(DependSuffix): ../../LidDrivenCavitySolver-master/src/LDCprogram_options.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCprogram_options.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCprogram_options.cpp$(DependSuffix) -MM ../../LidDrivenCavitySolver-master/src/LDCprogram_options.cpp

$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCprogram_options.cpp$(PreprocessSuffix): ../../LidDrivenCavitySolver-master/src/LDCprogram_options.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LDCprogram_options.cpp$(PreprocessSuffix) ../../LidDrivenCavitySolver-master/src/LDCprogram_options.cpp

$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavity.cpp$(ObjectSuffix): ../../LidDrivenCavitySolver-master/src/LidDrivenCavity.cpp $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavity.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/li/Desktop/LidDrivenCavitySolver-master/src/LidDrivenCavity.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavity.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavity.cpp$(DependSuffix): ../../LidDrivenCavitySolver-master/src/LidDrivenCavity.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavity.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavity.cpp$(DependSuffix) -MM ../../LidDrivenCavitySolver-master/src/LidDrivenCavity.cpp

$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavity.cpp$(PreprocessSuffix): ../../LidDrivenCavitySolver-master/src/LidDrivenCavity.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavity.cpp$(PreprocessSuffix) ../../LidDrivenCavitySolver-master/src/LidDrivenCavity.cpp

$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavitySolver.cpp$(ObjectSuffix): ../../LidDrivenCavitySolver-master/src/LidDrivenCavitySolver.cpp $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavitySolver.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/li/Desktop/LidDrivenCavitySolver-master/src/LidDrivenCavitySolver.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavitySolver.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavitySolver.cpp$(DependSuffix): ../../LidDrivenCavitySolver-master/src/LidDrivenCavitySolver.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavitySolver.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavitySolver.cpp$(DependSuffix) -MM ../../LidDrivenCavitySolver-master/src/LidDrivenCavitySolver.cpp

$(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavitySolver.cpp$(PreprocessSuffix): ../../LidDrivenCavitySolver-master/src/LidDrivenCavitySolver.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/up_up_LidDrivenCavitySolver-master_src_LidDrivenCavitySolver.cpp$(PreprocessSuffix) ../../LidDrivenCavitySolver-master/src/LidDrivenCavitySolver.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


