TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -std=c++11
QMAKE_LFLAGS += -static

SOURCES += \
        main.cpp \
    ldstat.cpp \
    cmdline.cpp \
    util.cpp \
    vcf.cpp

HEADERS += \
    cmdline.h \
    split.h \
    util.h \
    vcf.h
