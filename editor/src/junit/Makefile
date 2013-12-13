BASE = ..

include ${BASE}/Makefile.config

FILES = $(wildcard *.java)

CLASSES = ${FILES:.java=.class}

SUBDIRS = awtui extensions framework runner samples swingui tests textui

all: ${CLASSES} 

local-clean:
