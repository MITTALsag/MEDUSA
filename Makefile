# Variables
CXX = g++
CXXFLAGS = -g -std=c++17 -I. -DUSING_LzLog -DNO_GUI -DOS_LINUX -pthread
INCLUDES = -I./Eigen_3.3.7

# Liste des fichiers source de test
TEST = LzAsync/DispatchQ_UnitTest.cpp LzAsync/DispatchQ_Test.cpp

# Liste des fichiers source de LzServices (à ajouter)
LZSERVICES = LzServices/LzLog.cpp
DISPATCH = LzAsync/DispatchQ.cpp
DATA_THREAD = LzAsync/DataThread/DataThread.cpp
JOBS = LzAsync/Jobs/Jobs.cpp
GRAPH = LzAsync/GraphCycleThreads/GraphCycleThreads.cpp
LZTRIMODEL = $(wildcard LzTriModel/*.cpp)
MATH = $(wildcard LzMath/*.cpp)
GEOM = $(wildcard LzGeom/*.cpp)

# Liste des fichiers objets correspondants
OBJS = $(TEST:.cpp=.o) $(LZSERVICES:.cpp=.o) $(DISPATCH:.cpp=.o) $(JOBS:.cpp=.o) $(GRAPH:.cpp=.o) $(DATA_THREAD:.cpp=.o) $(LZTRIMODEL:.cpp=.o) $(MATH:.cpp=.o) $(GEOM:.cpp=.o)

# Nom de l'exécutable
EXEC = LzAsync/DispatchQ_Test

# Règle par défaut
all: $(EXEC)

# Création de l'exécutable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -I./Eigen_3.3.7 $^ -o $@

# Pattern rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Nettoyage des fichiers objets et exécutables
clean:
	rm -f $(OBJS) $(EXEC)


re : clean all
