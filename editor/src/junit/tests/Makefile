BASE = ../..

include ${BASE}/Makefile.config

FILES = $(wildcard *.java)

CLASSES = ${FILES:.java=.class}

SUBDIRS = extensions framework runner

all: ${CLASSES} 

local-clean:
