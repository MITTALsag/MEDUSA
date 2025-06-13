@INCLUPATH += LzTriModel # Creates the folder in project tree

HEADERS += $${CORE_SRC_DIR}/LzTriModel/LzLib_TriModel.h \
           $${CORE_SRC_DIR}/LzTriModel/Mesh.h \
           $${CORE_SRC_DIR}/LzTriModel/MeshTools.h \
           $${CORE_SRC_DIR}/LzTriModel/Triangle.h \
           $${CORE_SRC_DIR}/LzTriModel/MeshTopology.h \
           $${CORE_SRC_DIR}/LzTriModel/PatchComputer.h \
           $${CORE_SRC_DIR}/LzTriModel/FastDistComp.h \
           $${CORE_SRC_DIR}/LzTriModel/FastDistComp_Patching.h \
           $${CORE_SRC_DIR}/LzTriModel/Outliner.h \
           $${CORE_SRC_DIR}/LzTriModel/OutlinerTools.h \
           $${CORE_SRC_DIR}/LzTriModel/OutlinerTools_Meshing.h \
           $${CORE_SRC_DIR}/LzTriModel/TriangleMesher.h \
           $${CORE_SRC_DIR}/LzTriModel/DeformableGrid.h \
           $${CORE_SRC_DIR}/LzTriModel/MeshDecimator.h \

SOURCES += $${CORE_SRC_DIR}/LzTriModel/Mesh.cpp \
           $${CORE_SRC_DIR}/LzTriModel/MeshTools.cpp \
           $${CORE_SRC_DIR}/LzTriModel/MeshTools_Eigen.cpp \
           $${CORE_SRC_DIR}/LzTriModel/Triangle.cpp \
           $${CORE_SRC_DIR}/LzTriModel/MeshTopology.cpp \
           $${CORE_SRC_DIR}/LzTriModel/PatchComputer.cpp \
           $${CORE_SRC_DIR}/LzTriModel/FastDistComp.cpp \
           $${CORE_SRC_DIR}/LzTriModel/Outliner.cpp \
           $${CORE_SRC_DIR}/LzTriModel/OutlinerTools.cpp \
           $${CORE_SRC_DIR}/LzTriModel/OutlinerTools_Meshing.cpp \
           $${CORE_SRC_DIR}/LzTriModel/TriangleMesher.cpp \
           $${CORE_SRC_DIR}/LzTriModel/DeformableGrid.cpp \
           $${CORE_SRC_DIR}/LzTriModel/MeshDecimator.cpp \


