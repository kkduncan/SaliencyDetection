/*
 *  remviewer.cpp
 *  REM
 *
 *  Created by Kester Duncan on 12/26/10.
 *  Copyright 2010 Kester Duncan. All rights reserved.
 *
 */
#include <QtGui>

#include "remviewer.h"

REMView::REMView()
{
	imageLabel = new QLabel;
	imageLabel->setBackgroundRole(QPalette::Base);
	imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	imageLabel->setScaledContents(true);
	imageLabel->resize(500, 350);
	
	resultLabel = new QLabel;
	resultLabel->setBackgroundRole(QPalette::Base);
	resultLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	resultLabel->setScaledContents(true);
	resultLabel->resize(500, 350);
	
	mainLabel = new QLabel;
	mainLabel->setBackgroundRole(QPalette::Base);
	mainLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	mainLabel->setScaledContents(true);
	mainLabel->resize(1024, 728);
	
	scrollArea = new QScrollArea;
	scrollArea->setBackgroundRole(QPalette::Dark);
	scrollArea->setWidget(imageLabel);
	
	resultScrollArea = new QScrollArea;
	resultScrollArea->setBackgroundRole(QPalette::Dark);
	resultScrollArea->setWidget(resultLabel);
	
	QHBoxLayout *layout = new QHBoxLayout;
	layout->addWidget(scrollArea);
	layout->addWidget(resultScrollArea);
	mainLabel->setLayout(layout);
	
	setCentralWidget(mainLabel);	
	
	createActions();
	createMenus();
	
	QToolBar *fileToolBar = addToolBar(tr("File"));
	fileToolBar->addAction(openAct);
	
	QString iconName("remicon.png");
	QIcon icon(iconName);
	setWindowIcon(icon);
	setWindowTitle(tr("REM Saliency Viewer"));
	resize(1024, 728);
	
}

void REMView::open()
{
	fileName = QFileDialog::getOpenFileName(this,
											tr("Open File"), QDir::currentPath());
	if (!fileName.isEmpty()) {
		image = QImage(fileName);
		if (image.isNull()) {
			QMessageBox::information(this, tr("Image Viewer"),
									 tr("Cannot load %1. Accepted Image Formats - pgm, ppm, jpeg, png.").arg(fileName));
			return;
		}
		imageLabel->setPixmap(QPixmap::fromImage(image));
		scaleFactor = 1.0;
		
		printAct->setEnabled(true);
		fitToWindowAct->setEnabled(true);
		updateActions();
		
		if (!fitToWindowAct->isChecked())
			imageLabel->adjustSize();
	}
}

void REMView::print()
{
	Q_ASSERT(imageLabel->pixmap());
	QPrintDialog dialog(&printer, this);
	if (dialog.exec()) {
		QPainter painter(&printer);
		QRect rect = painter.viewport();
		QSize size = imageLabel->pixmap()->size();
		size.scale(rect.size(), Qt::KeepAspectRatio);
		painter.setViewport(rect.x(), rect.y(), size.width(), size.height());
		painter.setWindow(imageLabel->pixmap()->rect());
		painter.drawPixmap(0, 0, *imageLabel->pixmap());
	}
}

void REMView::computeSaliency() {
	//FILE *fp;
	//char *command;
	//sprintf(command, "/Users/kesterduncan/Documents/workspace/REM/bin/ComputeIterativeSaliency %s 25", fileName);
	//system(command);
	
	
	
}
void REMView::zoomIn()
{
	scaleImage(1.25);
}

void REMView::zoomOut()
{
	scaleImage(0.8);
}

void REMView::normalSize()
{
	imageLabel->adjustSize();
	scaleFactor = 1.0;
}

void REMView::fitToWindow()
{
	bool fitToWindow = fitToWindowAct->isChecked();
	scrollArea->setWidgetResizable(fitToWindow);
	if (!fitToWindow) {
		normalSize();
	}
	updateActions();
}

void REMView::about()
{
	QMessageBox::about(this, tr("About REM"),
					   tr("<p>The <b>Relational Entropy-Based Measure of Saliency</b> is a bottom-up" 
						  " measure of saliency based on the gradient direction and " 
						  " the Euclidean distance relationships between pixels.</p>"));
	
}

void REMView::createActions()
{
	openAct = new QAction(tr("&Open..."), this);
	openAct->setShortcut(tr("Ctrl+O"));
	connect(openAct, SIGNAL(triggered()), this, SLOT(open()));
	
	printAct = new QAction(tr("&Print..."), this);
	printAct->setShortcut(tr("Ctrl+P"));
	printAct->setEnabled(false);
	connect(printAct, SIGNAL(triggered()), this, SLOT(print()));
	
	exitAct = new QAction(tr("E&xit"), this);
	exitAct->setShortcut(tr("Ctrl+Q"));
	connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));
	
	zoomInAct = new QAction(tr("Zoom &In (25%)"), this);
	zoomInAct->setShortcut(tr("Ctrl++"));
	zoomInAct->setEnabled(false);
	connect(zoomInAct, SIGNAL(triggered()), this, SLOT(zoomIn()));
	
	zoomOutAct = new QAction(tr("Zoom &Out (25%)"), this);
	zoomOutAct->setShortcut(tr("Ctrl+-"));
	zoomOutAct->setEnabled(false);
	connect(zoomOutAct, SIGNAL(triggered()), this, SLOT(zoomOut()));
	
	normalSizeAct = new QAction(tr("&Normal Size"), this);
	normalSizeAct->setShortcut(tr("Ctrl+S"));
	normalSizeAct->setEnabled(false);
	connect(normalSizeAct, SIGNAL(triggered()), this, SLOT(normalSize()));
	
	fitToWindowAct = new QAction(tr("&Fit to Window"), this);
	fitToWindowAct->setEnabled(false);
	fitToWindowAct->setCheckable(true);
	fitToWindowAct->setShortcut(tr("Ctrl+F"));
	connect(fitToWindowAct, SIGNAL(triggered()), this, SLOT(fitToWindow()));
	
	aboutAct = new QAction(tr("&About"), this);
	connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));
	
}

void REMView::createMenus()
{
	fileMenu = new QMenu(tr("&File"), this);
	fileMenu->addAction(openAct);
	fileMenu->addAction(printAct);
	fileMenu->addSeparator();
	fileMenu->addAction(exitAct);
	
	viewMenu = new QMenu(tr("&View"), this);
	viewMenu->addAction(zoomInAct);
	viewMenu->addAction(zoomOutAct);
	viewMenu->addAction(normalSizeAct);
	viewMenu->addSeparator();
	viewMenu->addAction(fitToWindowAct);
	
	helpMenu = new QMenu(tr("&Help"), this);
	helpMenu->addAction(aboutAct);
	
	menuBar()->addMenu(fileMenu);
	menuBar()->addMenu(viewMenu);
	menuBar()->addMenu(helpMenu);
}

void REMView::updateActions()
{
	zoomInAct->setEnabled(!fitToWindowAct->isChecked());
	zoomOutAct->setEnabled(!fitToWindowAct->isChecked());
	normalSizeAct->setEnabled(!fitToWindowAct->isChecked());
}

void REMView::scaleImage(double factor)
{
	Q_ASSERT(imageLabel->pixmap());
	scaleFactor *= factor;
	imageLabel->resize(scaleFactor * imageLabel->pixmap()->size());
	
	adjustScrollBar(scrollArea->horizontalScrollBar(), factor);
	adjustScrollBar(scrollArea->verticalScrollBar(), factor);
	
	zoomInAct->setEnabled(scaleFactor < 3.0);
	zoomOutAct->setEnabled(scaleFactor > 0.333);
}

void REMView::adjustScrollBar(QScrollBar *scrollBar, double factor)
{
	scrollBar->setValue(int(factor * scrollBar->value()
							+ ((factor - 1) * scrollBar->pageStep()/2)));
}

