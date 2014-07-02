/*
 *  main.cpp
 *  REM
 *
 *  Created by Kester Duncan on 12/26/10.
 *  Copyright 2010 KDuncan. All rights reserved.
 *
 */
#include <QApplication>

#include "remviewer.h"

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);
	REMView remView;
	remView.show();
	return app.exec();
}