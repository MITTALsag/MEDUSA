@INCLUPATH += LzMath # Creates the folder in project tree

HEADERS += $${CORE_SRC_DIR}/LzMath/LzLib_Math.h \
           $${CORE_SRC_DIR}/LzMath/Matrix.h \
           $${CORE_SRC_DIR}/LzMath/Polynom.h \
           $${CORE_SRC_DIR}/LzMath/PolyMatrix.h \
           $${CORE_SRC_DIR}/LzMath/ToolBox.h \
           $${CORE_SRC_DIR}/LzMath/StringOperations.h \
           $${CORE_SRC_DIR}/LzMath/Graph.h \
           $${CORE_SRC_DIR}/LzMath/Graph_UnitTest.h \

SOURCES += $${CORE_SRC_DIR}/LzMath/Matrix.cpp \
           $${CORE_SRC_DIR}/LzMath/Polynom.cpp \
           $${CORE_SRC_DIR}/LzMath/PolyMatrix.cpp \
           $${CORE_SRC_DIR}/LzMath/ToolBox.cpp \
           $${CORE_SRC_DIR}/LzMath/ToolBox_Eigen.cpp \
           $${CORE_SRC_DIR}/LzMath/StringOperations.cpp \
           $${CORE_SRC_DIR}/LzMath/Graph.cpp \
           $${CORE_SRC_DIR}/LzMath/Graph_UnitTest.cpp \


