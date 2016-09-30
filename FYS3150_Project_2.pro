TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib

SOURCES += main.cpp \
    jacobi.cpp \
    interacting_case.cpp \
    non_interacting_case.cpp \
    unit_test.cpp

LIBS += -larmadillo -llapack -lblas

HEADERS += \
    jacobi.h \
    interacting_case.h \
    non_interacting_case.h \
    unit_test.h
