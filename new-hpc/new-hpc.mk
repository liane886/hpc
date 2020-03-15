##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=new-hpc
ConfigurationName      :=Debug
WorkspacePath          :=/home/li/Desktop/hpc
ProjectPath            :=/home/li/Desktop/hpc/new-hpc
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=li
Date                   :=14/03/20
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
ObjectsFileList        :="new-hpc.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := $(LibrarySwitch)blas $(LibrarySwitch)cblas $(LibrarySwitch)lapack 
ArLibs                 :=  "blas" "cblas" "lapack" 
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
Objects0=$(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(ObjectSuffix) $(IntermediateDirectory)/LidDrivenCavity.cpp$(ObjectSuffix) $(IntermediateDirectory)/Poisson.cpp$(ObjectSuffix) 



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
$(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(ObjectSuffix): LidDrivenCavitySolver.cpp $(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/li/Desktop/hpc/new-hpc/LidDrivenCavitySolver.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(DependSuffix): LidDrivenCavitySolver.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(DependSuffix) -MM LidDrivenCavitySolver.cpp

$(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(PreprocessSuffix): LidDrivenCavitySolver.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/LidDrivenCavitySolver.cpp$(PreprocessSuffix) LidDrivenCavitySolver.cpp

$(IntermediateDirectory)/LidDrivenCavity.cpp$(ObjectSuffix): LidDrivenCavity.cpp $(IntermediateDirectory)/LidDrivenCavity.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/li/Desktop/hpc/new-hpc/LidDrivenCavity.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/LidDrivenCavity.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/LidDrivenCavity.cpp$(DependSuffix): LidDrivenCavity.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/LidDrivenCavity.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/LidDrivenCavity.cpp$(DependSuffix) -MM LidDrivenCavity.cpp

$(IntermediateDirectory)/LidDrivenCavity.cpp$(PreprocessSuffix): LidDrivenCavity.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/LidDrivenCavity.cpp$(PreprocessSuffix) LidDrivenCavity.cpp

$(IntermediateDirectory)/Poisson.cpp$(ObjectSuffix): Poisson.cpp $(IntermediateDirectory)/Poisson.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/li/Desktop/hpc/new-hpc/Poisson.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Poisson.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Poisson.cpp$(DependSuffix): Poisson.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/Poisson.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/Poisson.cpp$(DependSuffix) -MM Poisson.cpp

$(IntermediateDirectory)/Poisson.cpp$(PreprocessSuffix): Poisson.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Poisson.cpp$(PreprocessSuffix) Poisson.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


