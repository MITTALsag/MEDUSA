@INCLUPATH += LzGeom # Creates the folder in project tree

HEADERS += $${CORE_SRC_DIR}/LzGeom/LzLib_Geom.h \
           $${CORE_SRC_DIR}/LzGeom/Coord3D.h \
           $${CORE_SRC_DIR}/LzGeom/Point3D.h \
           $${CORE_SRC_DIR}/LzGeom/Vector3D.h \
           $${CORE_SRC_DIR}/LzGeom/Line3D.h \
           $${CORE_SRC_DIR}/LzGeom/Plane3D.h \
           $${CORE_SRC_DIR}/LzGeom/RigidTr3D.h \
           $${CORE_SRC_DIR}/LzGeom/BBox.h \
           $${CORE_SRC_DIR}/LzGeom/Tree3D.h \

SOURCES += $${CORE_SRC_DIR}/LzGeom/Coord3D.cpp \
           $${CORE_SRC_DIR}/LzGeom/Point3D.cpp \
           $${CORE_SRC_DIR}/LzGeom/Vector3D.cpp \
           $${CORE_SRC_DIR}/LzGeom/Line3D.cpp \
           $${CORE_SRC_DIR}/LzGeom/Plane3D.cpp \
           $${CORE_SRC_DIR}/LzGeom/RigidTr3D.cpp \
           $${CORE_SRC_DIR}/LzGeom/BBox.cpp \
           $${CORE_SRC_DIR}/LzGeom/Tree3D.cpp \

