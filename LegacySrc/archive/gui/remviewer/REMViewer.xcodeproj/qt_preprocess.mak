#############################################################################
# Makefile for building: REMViewer.app/Contents/MacOS/REMViewer
# Generated by qmake (2.01a) (Qt 4.4.3) on: Sun Dec 26 12:37:42 2010
# Project:  REMViewer.pro
# Template: app
# Command: /usr/bin/qmake -spec /usr/local/Qt4.4/mkspecs/macx-xcode -macx -o REMViewer.xcodeproj/project.pbxproj REMViewer.pro
#############################################################################

MOC       = /Developer/Tools/Qt/moc
UIC       = /Developer/Tools/Qt/uic
LEX       = flex
LEXFLAGS  = 
YACC      = yacc
YACCFLAGS = -d
DEFINES       = -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED
INCPATH       = -I/usr/local/Qt4.4/mkspecs/macx-xcode -I. -I/Library/Frameworks/QtCore.framework/Versions/4/Headers -I/usr/include/QtCore -I/Library/Frameworks/QtGui.framework/Versions/4/Headers -I/usr/include/QtGui -I/usr/include -I. -I. -I. -I/usr/local/include -I/System/Library/Frameworks/CarbonCore.framework/Headers
DEL_FILE  = rm -f
MOVE      = mv -f

IMAGES = 
PARSERS =
preprocess: $(PARSERS) compilers
clean preprocess_clean: parser_clean compiler_clean

parser_clean:
mocclean: compiler_moc_header_clean compiler_moc_source_clean

mocables: compiler_moc_header_make_all compiler_moc_source_make_all

compilers: ./moc_remviewer.cpp
compiler_objective_c_make_all:
compiler_objective_c_clean:
compiler_moc_header_make_all: moc_remviewer.cpp
compiler_moc_header_clean:
	-$(DEL_FILE) moc_remviewer.cpp
moc_remviewer.cpp: remviewer.h
	/Developer/Tools/Qt/moc $(DEFINES) $(INCPATH) -D__APPLE__ -D__GNUC__ remviewer.h -o moc_remviewer.cpp

compiler_rcc_make_all:
compiler_rcc_clean:
compiler_image_collection_make_all: qmake_image_collection.cpp
compiler_image_collection_clean:
	-$(DEL_FILE) qmake_image_collection.cpp
compiler_moc_source_make_all:
compiler_moc_source_clean:
compiler_rez_source_make_all:
compiler_rez_source_clean:
compiler_uic_make_all:
compiler_uic_clean:
compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: compiler_moc_header_clean 

