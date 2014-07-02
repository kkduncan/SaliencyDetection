/****************************************************************************
** Meta object code from reading C++ file 'remviewer.h'
**
** Created: Tue Dec 28 19:32:50 2010
**      by: The Qt Meta Object Compiler version 59 (Qt 4.4.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "remviewer.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'remviewer.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.4.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_REMView[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
       9,    8,    8,    8, 0x08,
      16,    8,    8,    8, 0x08,
      24,    8,    8,    8, 0x08,
      33,    8,    8,    8, 0x08,
      43,    8,    8,    8, 0x08,
      56,    8,    8,    8, 0x08,
      70,    8,    8,    8, 0x08,
      78,    8,    8,    8, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_REMView[] = {
    "REMView\0\0open()\0print()\0zoomIn()\0"
    "zoomOut()\0normalSize()\0fitToWindow()\0"
    "about()\0computeSaliency()\0"
};

const QMetaObject REMView::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_REMView,
      qt_meta_data_REMView, 0 }
};

const QMetaObject *REMView::metaObject() const
{
    return &staticMetaObject;
}

void *REMView::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_REMView))
        return static_cast<void*>(const_cast< REMView*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int REMView::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: open(); break;
        case 1: print(); break;
        case 2: zoomIn(); break;
        case 3: zoomOut(); break;
        case 4: normalSize(); break;
        case 5: fitToWindow(); break;
        case 6: about(); break;
        case 7: computeSaliency(); break;
        }
        _id -= 8;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
