#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileDialog>
#include <QString>
#include <QFileInfo>
#include <QDir>
#include <QGridLayout>
#include <QTextEdit>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QLogValueAxis>
#include <QtCharts/QValueAxis>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QChart>
#include <QApplication>
#include <QProcess>
#include <QProgressBar>
#include <QString>

#include <ctime>
#include <fstream>
#include <vector>

#include "SAScore.h"
#include "SAXSParser.h"
#include "pdbsettings.h"
#include "about.h"

QT_CHARTS_USE_NAMESPACE
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    QChart *chart;
    QGridLayout* gridLayout;
    QGridLayout* gridLayout2D;

    double mult(std::vector<double>& data1, std::vector<double>& data2);
    void Rotate(std::vector<double>& data, double a);

    ~MainWindow();

private slots:
    void on_PDBchooseButon_clicked();

    void on_PipedCheckBox_clicked();

    void on_PipedChooseButon_clicked();

    void on_actionExit_triggered();

    void on_actionSettings_triggered();

    void on_ExitButton_clicked();

    void on_StartButton_clicked();

    void on_PlotButton_clicked();

    void on_actionAbout_triggered();

    void on_ExperCheckBox_clicked();

    void on_ExperChooseButon_clicked();

    void on_Plot2DButton_clicked();

    void on_PositionSlider_valueChanged(int value);

    void on_NormCoordEdit_editingFinished();

    void on_ButtonDown_clicked();

    void on_ButtonUp_clicked();

private:
    Ui::MainWindow *ui;
    QString ProjectPath_;
    QString PipedFilePath_;
    PDBSettings *Dialog;
    About *about;
    void Parsing();
    QDir Dir = QDir::currentPath();
    QString OutputFilePath_;
    void Log_state_(QString Process, QTextEdit *log, QProgressBar *Bar, int percent_state);
    void Output_File_Write_(QTextEdit *log);
    static bool XYZ_(const std::vector <int>& vec1, const std::vector <int>& vec2);
    static bool YZX_(const std::vector <int>& vec1, const std::vector <int>& vec2);
    static bool ZXY_(const std::vector <int>& vec1, const std::vector <int>& vec2);

    std::vector< std::vector<double> > Ie_;

    //Number of component of cube
    enum class CubeComponents_ {
        kX = 0,
        kY = 1,
        kZ = 2,
        kState = 3,
    };
};

#endif // MAINWINDOW_H
