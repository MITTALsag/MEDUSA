@INCLUPATH += LzServices # Creates the folder in project tree

HEADERS += $${CORE_SRC_DIR}/LzServices/LzLib_Services.h \
           $${CORE_SRC_DIR}/LzServices/LzLog.h \
           $${CORE_SRC_DIR}/LzServices/List.h \
           $${CORE_SRC_DIR}/LzServices/List_UnitTest.h \
           $${CORE_SRC_DIR}/LzServices/Array.h \
           $${CORE_SRC_DIR}/LzServices/Array_UnitTest.h \
           $${CORE_SRC_DIR}/LzServices/Array.inl \
           $${CORE_SRC_DIR}/LzServices/Vector.h \
           $${CORE_SRC_DIR}/LzServices/Vector.inl \
           $${CORE_SRC_DIR}/LzServices/HashTable.h \
           $${CORE_SRC_DIR}/LzServices/HashTable.inl \

SOURCES += $${CORE_SRC_DIR}/LzServices/LzLog.cpp \
           $${CORE_SRC_DIR}/LzServices/List_UnitTest.cpp \
           $${CORE_SRC_DIR}/LzServices/Array_UnitTest.cpp \


