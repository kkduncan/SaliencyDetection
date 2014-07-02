/*
 *  remviewer.h
 *  REM
 *
 *  Created by Kester Duncan on 12/26/10.
 *  Copyright 2010 KDuncan. All rights reserved.
 *
 */
#ifndef REMVIEWER_H
#define REMVIEWER_H

#include <QMainWindow>
#include <QHBoxLayout>
#include <QImage>
#include <QPrinter>
#include <cstdlib>
#include <cstdio>

class QAction;
class QLabel;
class QMenu;
class QImage;
class QScrollArea;
class QScrollBar;

class REMView : public QMainWindow
{
	Q_OBJECT
	
public:
	REMView();
	
	private slots:
	void open();
	void print();
	void zoomIn();
	void zoomOut();
	void normalSize();
	void fitToWindow();
	void about();
	void computeSaliency();
	
private:
	void createActions();
	void createMenus();
	void updateActions();
	void scaleImage(double factor);
	void adjustScrollBar(QScrollBar *scrollBar, double factor);
	
	QString fileName;
	QImage image;
	QLabel *mainLabel;
	QLabel *imageLabel;
	QLabel *resultLabel;
	QScrollArea *scrollArea;
	QScrollArea *resultScrollArea;
	double scaleFactor;
	
	QPrinter printer;
	
	QAction *openAct;
	QAction *printAct;
	QAction *exitAct;
	QAction *zoomInAct;
	QAction *zoomOutAct;
	QAction *normalSizeAct;
	QAction *fitToWindowAct;
	QAction *aboutAct;
	//QAction *computeSaliencyAct;
	
	QMenu *fileMenu;
	QMenu *viewMenu;
	QMenu *helpMenu;
};

#endif

